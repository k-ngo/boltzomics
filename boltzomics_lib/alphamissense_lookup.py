"""
AlphaMissense lookup and prioritization helpers for mutation discovery.

This module is intentionally standalone so the integration can be removed by:
1) Removing the import/call sites in UI code, and
2) Deleting this file.
"""

from __future__ import annotations

import gzip
import os
from functools import lru_cache
from typing import Any, Dict, List, Optional, Sequence, Tuple


DEFAULT_AM_GZ = os.path.join("alphamissense", "AlphaMissense_aa_substitutions.by_uniprot.tsv.gz")
DEFAULT_AM_TBI = DEFAULT_AM_GZ + ".tbi"


def _severity_to_score(severity: Optional[str]) -> float:
    value = str(severity or "").strip().lower()
    if value == "high":
        return 1.0
    if value == "moderate":
        return 0.6
    if value == "low":
        return 0.25
    return 0.0


def _priority_tier(score: float) -> str:
    if score >= 0.75:
        return "HIGH"
    if score >= 0.5:
        return "MEDIUM"
    return "LOW"


def _normalize_variant_code(code: Optional[str]) -> str:
    return str(code or "").strip().upper()


def _resolve_canonical_variant_code(entry: Dict[str, Any]) -> str:
    canonical_code = _normalize_variant_code(entry.get("canonical_code"))
    if canonical_code:
        return canonical_code

    wt = _normalize_variant_code(entry.get("wt_residue"))
    mut = _normalize_variant_code(entry.get("mut_residue"))
    try:
        canonical_pos = int(entry.get("canonical_position"))
    except Exception:
        canonical_pos = None
    if wt and mut and canonical_pos and canonical_pos > 0:
        return f"{wt[0]}{int(canonical_pos)}{mut[0]}"

    return _normalize_variant_code(entry.get("mutation_code"))


class AlphaMissenseLookup:
    """
    Fast indexed AlphaMissense reader using tabix index.

    Expected columns:
    #uniprot_id, position, protein_variant, am_pathogenicity, am_class
    """

    def __init__(self, tsv_gz_path: str = DEFAULT_AM_GZ, tbi_path: Optional[str] = None):
        self.tsv_gz_path = tsv_gz_path
        self.tbi_path = tbi_path or (tsv_gz_path + ".tbi")
        self._tabix = None
        self._lookup_backend = "none"
        self._cache: Dict[str, Dict[str, Dict[str, Any]]] = {}
        self._available = False
        self._status_reason = "not_initialized"
        self._initialize()

    def _initialize(self) -> None:
        if not os.path.exists(self.tsv_gz_path):
            self._status_reason = f"missing_file:{self.tsv_gz_path}"
            return

        has_index = os.path.exists(self.tbi_path)
        try:
            import pysam  # type: ignore
            has_pysam = True
        except Exception as exc:
            has_pysam = False
            pysam_error = str(exc)

        if has_pysam and has_index:
            try:
                self._tabix = pysam.TabixFile(self.tsv_gz_path)
                self._lookup_backend = "tabix"
                self._available = True
                self._status_reason = "ok"
                return
            except Exception as exc:
                self._status_reason = f"tabix_open_failed:{exc}"

        # Fallback path: scan gzipped table and cache by UniProt ID.
        self._lookup_backend = "gzip_scan"
        self._available = True
        if not has_pysam:
            self._status_reason = f"fallback_gzip_scan_no_pysam ({pysam_error})"
        elif not has_index:
            self._status_reason = f"fallback_gzip_scan_missing_index:{self.tbi_path}"
        else:
            self._status_reason = "fallback_gzip_scan"

    @property
    def available(self) -> bool:
        return self._available

    @property
    def status_reason(self) -> str:
        return self._status_reason

    def _candidate_uniprot_ids(self, uniprot_id: str) -> List[str]:
        raw = str(uniprot_id or "").strip()
        if not raw:
            return []
        candidates = [raw]
        if "-" in raw:
            base = raw.split("-", 1)[0]
            if base:
                candidates.append(base)
        # Keep order, remove duplicates.
        unique: List[str] = []
        seen = set()
        for item in candidates:
            if item not in seen:
                unique.append(item)
                seen.add(item)
        return unique

    def _fetch_uniprot_rows(self, uniprot_id: str) -> Dict[str, Dict[str, Any]]:
        if not self.available:
            return {}
        if self._lookup_backend == "tabix":
            return self._fetch_uniprot_rows_tabix(uniprot_id)
        return self._fetch_uniprot_rows_gzip(uniprot_id)

    def _parse_line_to_record(self, parts: List[str]) -> Optional[Tuple[str, Dict[str, Any]]]:
        if len(parts) < 5:
            return None
        variant = _normalize_variant_code(parts[2])
        if not variant:
            return None
        try:
            position = int(parts[1])
        except Exception:
            position = None
        try:
            am_score = float(parts[3])
        except Exception:
            am_score = None
        am_class = str(parts[4]).strip().lower() if parts[4] else None
        return (
            variant,
            {
                "uniprot_id": str(parts[0]).strip(),
                "position": position,
                "protein_variant": variant,
                "am_pathogenicity": am_score,
                "am_class": am_class,
            },
        )

    def _fetch_uniprot_rows_tabix(self, uniprot_id: str) -> Dict[str, Dict[str, Any]]:
        if self._tabix is None:
            return {}

        try:
            lines = self._tabix.fetch(uniprot_id)
        except TypeError:
            try:
                lines = self._tabix.fetch(reference=uniprot_id)
            except Exception:
                return {}
        except Exception:
            return {}

        records: Dict[str, Dict[str, Any]] = {}
        for line in lines:
            if not line or line.startswith("#"):
                continue
            parsed = self._parse_line_to_record(line.rstrip("\n").split("\t"))
            if parsed is None:
                continue
            variant, record = parsed
            records[variant] = record
        return records

    def _fetch_uniprot_rows_gzip(self, uniprot_id: str) -> Dict[str, Dict[str, Any]]:
        matches = self._fetch_uniprot_rows_gzip_bulk([uniprot_id])
        return matches.get(uniprot_id, {})

    def _fetch_uniprot_rows_gzip_bulk(self, uniprot_ids: Sequence[str]) -> Dict[str, Dict[str, Dict[str, Any]]]:
        targets = {str(item).strip() for item in uniprot_ids if str(item).strip()}
        out: Dict[str, Dict[str, Dict[str, Any]]] = {item: {} for item in targets}
        if not targets:
            return out

        try:
            with gzip.open(self.tsv_gz_path, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    uid = str(parts[0]).strip()
                    if uid not in targets:
                        continue
                    parsed = self._parse_line_to_record(parts)
                    if parsed is None:
                        continue
                    variant, record = parsed
                    out[uid][variant] = record
        except Exception:
            return out
        return out

    def get_variant_map(self, uniprot_id: str) -> Dict[str, Dict[str, Any]]:
        candidates = self._candidate_uniprot_ids(uniprot_id)
        if self._lookup_backend == "gzip_scan":
            missing = [candidate for candidate in candidates if candidate not in self._cache]
            if missing:
                fetched = self._fetch_uniprot_rows_gzip_bulk(missing)
                for candidate in missing:
                    self._cache[candidate] = fetched.get(candidate, {})
        for candidate in candidates:
            if candidate not in self._cache:
                self._cache[candidate] = self._fetch_uniprot_rows(candidate)
            cached = self._cache.get(candidate) or {}
            if cached:
                return cached
        return {}

    def lookup_variant(self, uniprot_id: str, protein_variant: str) -> Optional[Dict[str, Any]]:
        variant = _normalize_variant_code(protein_variant)
        if not variant:
            return None
        variant_map = self.get_variant_map(uniprot_id)
        return variant_map.get(variant)


@lru_cache(maxsize=1)
def get_default_alphamissense_lookup() -> AlphaMissenseLookup:
    return AlphaMissenseLookup(DEFAULT_AM_GZ, DEFAULT_AM_TBI)


def annotate_mutations_with_alphamissense(
    mutations: List[Dict[str, Any]],
    use_prioritization: bool = True,
    min_am_score: Optional[float] = None,
    allowed_am_classes: Optional[Sequence[str]] = None,
    top_n: Optional[int] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Add AlphaMissense annotations and optional ranking/filtering.

    Returns:
      (annotated_mutations, summary)
    """
    if not mutations:
        return [], {
            "available": False,
            "reason": "empty_input",
            "matched": 0,
            "total": 0,
        }

    lookup = get_default_alphamissense_lookup()
    summary = {
        "available": lookup.available,
        "reason": lookup.status_reason,
        "matched": 0,
        "total": len(mutations),
    }

    allowed_classes_set = None
    if allowed_am_classes:
        allowed_classes_set = {str(item).strip().lower() for item in allowed_am_classes if str(item).strip()}

    annotated: List[Dict[str, Any]] = []
    for row in mutations:
        entry = dict(row)
        canonical_code = _resolve_canonical_variant_code(entry)
        accession = str(entry.get("uniprot_accession") or "").strip()
        am_record = lookup.lookup_variant(accession, canonical_code) if lookup.available else None
        if am_record:
            summary["matched"] += 1
            entry["am_pathogenicity"] = am_record.get("am_pathogenicity")
            entry["am_class"] = am_record.get("am_class")
            entry["am_lookup_status"] = "matched"
        else:
            entry["am_pathogenicity"] = None
            entry["am_class"] = None
            entry["am_lookup_status"] = "not_found" if lookup.available else "unavailable"

        severity_score = _severity_to_score(entry.get("severity"))
        am_score = entry.get("am_pathogenicity")

        # --- SIFT contribution (0-1, inverted: lower SIFT = more damaging) ---
        sift_component = 0.0
        has_sift = False
        raw_sift = entry.get("sift_score")
        if raw_sift is not None:
            try:
                sift_val = float(raw_sift)
                sift_component = 1.0 - sift_val  # invert: SIFT 0 = damaging → 1.0
                has_sift = True
            except (ValueError, TypeError):
                pass

        # --- PolyPhen contribution (0-1, already higher = more damaging) ---
        polyphen_component = 0.0
        has_polyphen = False
        raw_polyphen = entry.get("polyphen_score")
        if raw_polyphen is not None:
            try:
                polyphen_component = float(raw_polyphen)
                has_polyphen = True
            except (ValueError, TypeError):
                pass

        # --- Clinical significance bonus ---
        clin_sig = str(entry.get("clinical_significance") or "").strip().lower()
        clin_bonus = 0.0
        if any(term in clin_sig for term in ("pathogenic", "likely pathogenic", "risk factor")):
            clin_bonus = 0.15
        elif "uncertain" in clin_sig:
            clin_bonus = 0.05

        # --- Build weighted priority score dynamically ---
        components: list = []  # (weight, value) pairs
        sources: list = []

        if am_score is not None:
            components.append((0.40, float(am_score)))
            sources.append("alphamissense")
        if has_sift:
            components.append((0.15, sift_component))
            sources.append("sift")
        if has_polyphen:
            components.append((0.15, polyphen_component))
            sources.append("polyphen")

        # Severity (legacy) always contributes
        remaining_weight = 1.0 - sum(w for w, _ in components)
        components.append((remaining_weight, severity_score))
        sources.append("legacy")

        priority_score = sum(w * v for w, v in components) + clin_bonus
        priority_score = max(0.0, min(1.0, float(priority_score)))
        priority_source = "+".join(sources)

        entry["priority_score"] = priority_score
        entry["priority_tier"] = _priority_tier(priority_score)
        entry["priority_source"] = priority_source
        annotated.append(entry)

    filtered = annotated
    if min_am_score is not None:
        filtered = [
            item for item in filtered
            if item.get("am_pathogenicity") is not None and float(item["am_pathogenicity"]) >= float(min_am_score)
        ]

    if allowed_classes_set is not None and len(allowed_classes_set) > 0:
        filtered = [
            item for item in filtered
            if str(item.get("am_class") or "").strip().lower() in allowed_classes_set
        ]

    if use_prioritization:
        filtered.sort(
            key=lambda item: (
                float(item.get("priority_score", 0.0)),
                float(item.get("evidence_count", 0.0)),
            ),
            reverse=True,
        )

    if top_n is not None and top_n > 0:
        filtered = filtered[: int(top_n)]

    summary["returned"] = len(filtered)
    return filtered, summary
