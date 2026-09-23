#!/usr/bin/env python3
"""
Example: MSA Caching for Accelerated Mutation Screening

This script demonstrates how the MSA caching system accelerates
pharmacogenomic screening of protein variants.

The key insight is that point mutations don't significantly change
the evolutionary context of a protein - the same homologous sequences
will be found during MSA generation. Therefore, we can:

1. Generate MSA once for the wild-type protein
2. Reuse it for all mutant variants
3. Achieve 95-99% reduction in MSA computation time

Run this example:
    python examples/msa_caching_example.py
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from msa_cache import MSACache, get_msa_cache


def demonstrate_msa_caching():
    """Demonstrate the MSA caching workflow."""

    # Example protein sequences (CYP3A4 fragment for demo)
    wt_sequence = "MALIPDLAMETWLLLAVSLVLLYLYGTHSHGLFKKLGIPGPTPLPFLGNILSYHKGF"
    k262r_mutant = "MALIPDLAMETWLLLAVSLVLLYLYGTHSHGLFKKLGIPGPTPLPFLGNILSYHRGF"  # K→R

    print("=" * 60)
    print("MSA Caching Demonstration")
    print("=" * 60)

    # Initialize cache
    cache = get_msa_cache(cache_dir="msa_cache_example")
    print(f"\nCache directory: {cache.cache_dir}")

    # Check cache statistics
    stats = cache.get_cache_stats()
    print(f"Cache entries: {stats['total_entries']}")
    print(f"Cache size: {stats['total_size_mb']} MB")

    # Compute sequence hashes
    wt_hash = MSACache.compute_sequence_hash(wt_sequence)
    mutant_hash = MSACache.compute_sequence_hash(k262r_mutant)

    print(f"\nWT sequence hash: {wt_hash[:16]}...")
    print(f"Mutant sequence hash: {mutant_hash[:16]}...")

    # Check for cached MSA
    print("\n--- Checking for cached MSA ---")

    wt_msa = cache.get_cached_msa(wt_sequence)
    if wt_msa:
        print(f"WT MSA found in cache: {wt_msa}")
    else:
        print("WT MSA not in cache - would need to generate")

    # For mutants, we can fall back to WT MSA
    mutant_msa = cache.get_cached_msa_for_mutant(wt_sequence, k262r_mutant)
    if mutant_msa:
        print(f"Mutant can use cached MSA: {mutant_msa}")
    else:
        print("No cached MSA available for mutant")

    # Demonstrate the workflow
    print("\n--- Mutation Screening Workflow ---")
    print("""
    Traditional Approach (without caching):
    ----------------------------------------
    For 10 mutants × 20 drugs = 200 predictions:
    - Each prediction generates its own MSA (~45 seconds)
    - Total MSA time: 200 × 45s = 2.5 hours

    With MSA Caching:
    -----------------
    1. First prediction (WT + first drug):
       - Generate MSA (~45 seconds)
       - Cache the MSA for reuse

    2. Subsequent predictions (all other combinations):
       - Reuse cached MSA (0 seconds MSA time)
       - Total MSA time: 1 × 45s = 45 seconds

    Speedup: ~97% reduction in MSA computation time!
    """)

    # Show how to integrate with Boltz prediction
    print("--- Integration with Boltz Prediction ---")
    print("""
    # In your screening code:
    from msa_cache import get_msa_cache

    cache = get_msa_cache()

    # Check if MSA is cached
    msa_path = cache.get_cached_msa_for_mutant(wt_seq, mutant_seq)

    if msa_path:
        # Use cached MSA - much faster!
        yaml_content = {
            "sequences": [{
                "protein": {
                    "id": "A",
                    "sequence": mutant_seq,
                    "msa": msa_path  # <-- Point to cached MSA
                }
            }, ...]
        }
        # Run Boltz WITHOUT --use_msa_server flag
        run_boltz_prediction(..., use_cached_msa=True)
    else:
        # Generate MSA normally, then cache it
        run_boltz_prediction(...)
        cache.cache_msa_from_boltz_output(wt_seq, boltz_output_dir)
    """)

    print("=" * 60)
    print("See NOVELTY.md for full documentation")
    print("=" * 60)


if __name__ == "__main__":
    demonstrate_msa_caching()
