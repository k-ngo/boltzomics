"""
MSA (Multiple Sequence Alignment) Cache Manager for BoltzOmics

This module implements a caching strategy for MSA files to dramatically accelerate
pharmacogenomic screening of protein variants. The key insight is that point mutations
do not significantly change the evolutionary context of a protein - the same homologous
sequences will be found during database searches.

Scientific Justification:
-------------------------
1. MSA captures evolutionary covariance between residue positions, not exact sequences
2. Point mutations don't change protein family membership - same homologs are found
3. Boltz-2 uses MSA for co-evolutionary signal which is conserved across variants
4. This is consistent with AlphaFold2/3 approaches where MSA provides family-level info

Performance Impact:
------------------
- MSA generation typically takes 30-60 seconds per query (network + compute)
- For mutation screening: N mutants × M drugs = N×M redundant MSA generations
- With caching: 1 MSA generation, reused across all combinations
- Typical savings: 95-99% reduction in MSA computation time

Usage:
------
    from msa_cache import MSACache

    cache = MSACache()

    # Check if MSA exists for this protein
    msa_path = cache.get_cached_msa(protein_sequence)

    if msa_path is None:
        # Run Boltz with MSA generation, then cache
        run_boltz_prediction(...)  # generates MSA
        cache.cache_msa_from_boltz_output(protein_sequence, boltz_output_dir)
    else:
        # Use cached MSA in YAML
        create_yaml_with_msa(protein_sequence, msa_path, ...)
"""

import os
import hashlib
import logging
import shutil
import json
import glob
from datetime import datetime
from typing import Optional, Dict, List, Tuple
import csv

logger = logging.getLogger(__name__)


class MSACache:
    """
    Manages MSA file caching for efficient variant screening.

    Cache Structure:
    ----------------
    msa_cache/
    ├── index.json              # Maps sequence hashes to metadata
    ├── <hash_prefix>/
    │   └── <full_hash>/
    │       ├── msa_paired.csv      # Paired MSA file
    │       ├── msa_unpaired.csv    # Unpaired MSA file (if exists)
    │       └── metadata.json       # Creation time, source, etc.
    """

    def __init__(self, cache_dir: str = "msa_cache"):
        """
        Initialize the MSA cache.

        Args:
            cache_dir: Root directory for MSA cache storage
        """
        self.cache_dir = cache_dir
        self.index_path = os.path.join(cache_dir, "index.json")
        self._ensure_cache_dir()
        self._index = self._load_index()

    def _ensure_cache_dir(self) -> None:
        """Create cache directory if it doesn't exist."""
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

    def _load_index(self) -> Dict:
        """Load the cache index from disk."""
        if os.path.exists(self.index_path):
            try:
                with open(self.index_path, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                return {"sequences": {}, "version": "1.0"}
        return {"sequences": {}, "version": "1.0"}

    def _save_index(self) -> None:
        """Save the cache index to disk."""
        with open(self.index_path, 'w') as f:
            json.dump(self._index, f, indent=2)

    @staticmethod
    def compute_sequence_hash(sequence: str) -> str:
        """
        Compute a deterministic hash for a protein sequence.

        The hash is computed on the uppercase, stripped sequence to ensure
        consistent lookups regardless of formatting differences.

        Args:
            sequence: Protein amino acid sequence

        Returns:
            SHA-256 hash of the normalized sequence
        """
        # Normalize: uppercase, remove whitespace, handle multi-chain
        normalized = sequence.upper().strip()
        # Remove any non-amino-acid characters except chain separator ':'
        normalized = ''.join(c for c in normalized if c.isalpha() or c == ':')
        return hashlib.sha256(normalized.encode('utf-8')).hexdigest()

    def _get_cache_path(self, seq_hash: str) -> str:
        """Get the cache directory path for a sequence hash."""
        # Use first 2 characters as prefix for better filesystem distribution
        prefix = seq_hash[:2]
        return os.path.join(self.cache_dir, prefix, seq_hash)

    def get_cached_msa(self, sequence: str) -> Optional[str]:
        """
        Check if a cached MSA exists for the given protein sequence.

        Args:
            sequence: Protein amino acid sequence

        Returns:
            Path to cached MSA file if exists, None otherwise
        """
        seq_hash = self.compute_sequence_hash(sequence)

        if seq_hash not in self._index.get("sequences", {}):
            return None

        cache_path = self._get_cache_path(seq_hash)
        msa_file = os.path.join(cache_path, "msa_paired.csv")

        if os.path.exists(msa_file):
            return msa_file

        # Check for .a3m format as fallback
        msa_a3m = os.path.join(cache_path, "msa_paired.a3m")
        if os.path.exists(msa_a3m):
            return msa_a3m

        # MSA file missing, remove from index
        del self._index["sequences"][seq_hash]
        self._save_index()
        return None

    def get_cached_msa_for_mutant(self, wt_sequence: str, mutant_sequence: str) -> Optional[str]:
        """
        Get cached MSA for a mutant, falling back to WT MSA if available.

        For point mutations, we can reuse the WT MSA since the evolutionary
        context (homologous sequences) is essentially identical.

        Args:
            wt_sequence: Wild-type protein sequence
            mutant_sequence: Mutant protein sequence

        Returns:
            Path to cached MSA file, preferring mutant-specific, falling back to WT
        """
        # First check for mutant-specific MSA
        mutant_msa = self.get_cached_msa(mutant_sequence)
        if mutant_msa:
            return mutant_msa

        # Fall back to WT MSA (scientifically valid for point mutations)
        wt_msa = self.get_cached_msa(wt_sequence)
        if wt_msa:
            logger.warning(
                "No mutant-specific MSA cached (hash=%s); falling back to WT MSA (hash=%s)",
                self.compute_sequence_hash(mutant_sequence)[:12],
                self.compute_sequence_hash(wt_sequence)[:12],
            )
        return wt_msa

    def cache_msa_from_boltz_output(
        self,
        sequence: str,
        boltz_output_dir: str,
        source_info: Optional[Dict] = None
    ) -> Optional[str]:
        """
        Cache MSA files from a Boltz prediction output directory.

        Boltz stores MSA files in: boltz_results_<name>/msa/<name>_0.csv

        Args:
            sequence: Protein sequence that was used for prediction
            boltz_output_dir: Path to boltz_results_* directory
            source_info: Optional metadata about the source prediction

        Returns:
            Path to cached MSA file, or None if MSA not found
        """
        seq_hash = self.compute_sequence_hash(sequence)

        # Find MSA files in Boltz output
        msa_dir = os.path.join(boltz_output_dir, "msa")
        if not os.path.exists(msa_dir):
            return None

        # Look for CSV MSA files (Boltz format)
        msa_files = glob.glob(os.path.join(msa_dir, "*_0.csv"))
        if not msa_files:
            # Try .a3m format
            msa_files = glob.glob(os.path.join(msa_dir, "*.a3m"))

        if not msa_files:
            return None

        # Create cache directory
        cache_path = self._get_cache_path(seq_hash)
        os.makedirs(cache_path, exist_ok=True)

        # Copy MSA file to cache
        source_msa = msa_files[0]
        if source_msa.endswith('.csv'):
            dest_msa = os.path.join(cache_path, "msa_paired.csv")
        else:
            dest_msa = os.path.join(cache_path, "msa_paired.a3m")

        shutil.copy2(source_msa, dest_msa)

        # Also copy unpaired MSA if exists
        unpaired_files = glob.glob(os.path.join(msa_dir, "*_1.csv"))
        if unpaired_files:
            shutil.copy2(unpaired_files[0], os.path.join(cache_path, "msa_unpaired.csv"))

        # Save metadata
        metadata = {
            "sequence_hash": seq_hash,
            "sequence_length": len(sequence.replace(":", "")),
            "cached_at": datetime.now().isoformat(),
            "source_dir": boltz_output_dir,
            "source_info": source_info or {}
        }
        with open(os.path.join(cache_path, "metadata.json"), 'w') as f:
            json.dump(metadata, f, indent=2)

        # Update index
        self._index["sequences"][seq_hash] = {
            "cached_at": metadata["cached_at"],
            "sequence_length": metadata["sequence_length"]
        }
        self._save_index()

        return dest_msa

    def create_mutant_msa(
        self,
        wt_sequence: str,
        mutant_sequence: str,
        mutations: List[Tuple[str, int, str, str]]  # [(chain, pos, wt_res, mut_res), ...]
    ) -> Optional[str]:
        """
        Create a mutant-specific MSA by modifying the query sequence in cached WT MSA.

        This is more accurate than directly reusing WT MSA as it ensures the
        query sequence in the MSA matches the actual mutant sequence.

        Args:
            wt_sequence: Wild-type protein sequence
            mutant_sequence: Mutant protein sequence
            mutations: List of mutations as (chain_id, position, wt_residue, mut_residue)

        Returns:
            Path to created mutant MSA, or None if WT MSA not cached
        """
        wt_msa_path = self.get_cached_msa(wt_sequence)
        if not wt_msa_path:
            return None

        mutant_hash = self.compute_sequence_hash(mutant_sequence)

        # Check if mutant MSA already exists
        existing = self.get_cached_msa(mutant_sequence)
        if existing:
            return existing

        # Create mutant cache directory
        cache_path = self._get_cache_path(mutant_hash)
        os.makedirs(cache_path, exist_ok=True)

        # Validate source MSA before modification
        with open(wt_msa_path, "r") as f:
            source_line_count = sum(1 for _ in f)
        if source_line_count < 3:
            raise ValueError(
                f"Source WT MSA appears corrupt: only {source_line_count} lines "
                f"in {wt_msa_path} (expected header + query + homologs)"
            )

        # Read WT MSA and modify query sequence
        if wt_msa_path.endswith('.csv'):
            mutant_msa_path = os.path.join(cache_path, "msa_paired.csv")
            self._modify_csv_msa_query(wt_msa_path, mutant_msa_path, wt_sequence, mutant_sequence)
        else:
            mutant_msa_path = os.path.join(cache_path, "msa_paired.a3m")
            self._modify_a3m_msa_query(wt_msa_path, mutant_msa_path, wt_sequence, mutant_sequence)

        # Validate output MSA has same line count as source
        with open(mutant_msa_path, "r") as f:
            output_line_count = sum(1 for _ in f)
        if output_line_count != source_line_count:
            logger.error(
                "Mutant MSA line count mismatch: source=%d, output=%d. Removing corrupt output.",
                source_line_count, output_line_count,
            )
            os.remove(mutant_msa_path)
            raise ValueError(
                f"Mutant MSA line count mismatch: source={source_line_count}, "
                f"output={output_line_count}"
            )

        # Save metadata
        metadata = {
            "sequence_hash": mutant_hash,
            "sequence_length": len(mutant_sequence.replace(":", "")),
            "cached_at": datetime.now().isoformat(),
            "derived_from_wt": self.compute_sequence_hash(wt_sequence),
            "mutations": [{"chain": m[0], "position": m[1], "wt": m[2], "mut": m[3]} for m in mutations]
        }
        with open(os.path.join(cache_path, "metadata.json"), 'w') as f:
            json.dump(metadata, f, indent=2)

        # Update index
        self._index["sequences"][mutant_hash] = {
            "cached_at": metadata["cached_at"],
            "sequence_length": metadata["sequence_length"],
            "derived_from": metadata["derived_from_wt"]
        }
        self._save_index()

        return mutant_msa_path

    def _modify_csv_msa_query(
        self,
        source_path: str,
        dest_path: str,
        wt_sequence: str,
        mutant_sequence: str
    ) -> None:
        """Modify the query sequence (row 0) in a CSV MSA file."""
        with open(source_path, 'r', newline='') as src, \
             open(dest_path, 'w', newline='') as dst:
            reader = csv.reader(src)
            writer = csv.writer(dst)

            # Copy header
            header = next(reader)
            writer.writerow(header)

            # Modify first data row (query sequence)
            first_row = next(reader)
            if len(first_row) >= 2:
                # Row format: key, sequence
                first_row[1] = mutant_sequence.replace(":", "")
            writer.writerow(first_row)

            # Copy remaining rows unchanged
            for row in reader:
                writer.writerow(row)

    def _modify_a3m_msa_query(
        self,
        source_path: str,
        dest_path: str,
        wt_sequence: str,
        mutant_sequence: str
    ) -> None:
        """Modify the query sequence in an A3M MSA file."""
        with open(source_path, 'r') as src, open(dest_path, 'w') as dst:
            lines = src.readlines()

            # A3M format: >header\nsequence\n...
            # First sequence is the query
            i = 0
            while i < len(lines):
                line = lines[i]
                if line.startswith('>'):
                    dst.write(line)
                    i += 1
                    # Next line is sequence
                    if i < len(lines) and not lines[i].startswith('>'):
                        if i == 1:  # First sequence (query)
                            dst.write(mutant_sequence.replace(":", "") + '\n')
                        else:
                            dst.write(lines[i])
                        i += 1
                else:
                    dst.write(line)
                    i += 1

    def get_cache_stats(self) -> Dict:
        """Get statistics about the MSA cache."""
        total_entries = len(self._index.get("sequences", {}))

        # Calculate total size
        total_size = 0
        for root, dirs, files in os.walk(self.cache_dir):
            for f in files:
                total_size += os.path.getsize(os.path.join(root, f))

        return {
            "total_entries": total_entries,
            "total_size_mb": round(total_size / (1024 * 1024), 2),
            "cache_dir": self.cache_dir
        }

    def clear_cache(self) -> None:
        """Clear all cached MSA files."""
        if os.path.exists(self.cache_dir):
            shutil.rmtree(self.cache_dir)
        self._ensure_cache_dir()
        self._index = {"sequences": {}, "version": "1.0"}
        self._save_index()


# Global cache instance
_msa_cache: Optional[MSACache] = None


def get_msa_cache(cache_dir: str = "msa_cache") -> MSACache:
    """Get or create the global MSA cache instance."""
    global _msa_cache
    if _msa_cache is None:
        _msa_cache = MSACache(cache_dir)
    return _msa_cache


def should_use_msa_cache(
    protein_sequence: str,
    wt_sequence: Optional[str] = None,
    enable_cache: bool = True
) -> Tuple[bool, Optional[str]]:
    """
    Determine if MSA cache should be used for a prediction.

    Args:
        protein_sequence: The protein sequence to predict
        wt_sequence: Optional WT sequence (for mutant predictions)
        enable_cache: Whether caching is enabled

    Returns:
        Tuple of (should_use_cache, msa_path_if_cached)
    """
    if not enable_cache:
        return False, None

    cache = get_msa_cache()

    # For mutants, try to use WT MSA
    if wt_sequence and wt_sequence != protein_sequence:
        msa_path = cache.get_cached_msa_for_mutant(wt_sequence, protein_sequence)
    else:
        msa_path = cache.get_cached_msa(protein_sequence)

    return msa_path is not None, msa_path
