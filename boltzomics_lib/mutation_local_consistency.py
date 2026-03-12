"""
Mutation-local consistency annotation for BoltzOmics.

This module is intentionally standalone so the feature can be enabled/disabled
without affecting core screening execution code.
"""

from __future__ import annotations

import os
import re
from typing import Any, Dict, List, Optional, Tuple

try:
    from structure_refinement import compare_wt_mutant_quick
except Exception:
    compare_wt_mutant_quick = None


def _is_wt_name(name: str) -> bool:
    upper = str(name).strip().upper()
    return upper in {"WT", "WILDTYPE", "WILD_TYPE"} or upper.endswith("_WT") or upper.endswith("-WT")


def _extract_mutation_positions(protein_name: str) -> List[int]:
    if not protein_name:
        return []
    matches = re.findall(r"[A-Za-z](\d+)[A-Za-z]", str(protein_name))
    positions: List[int] = []
    for token in matches:
        try:
            pos = int(token)
        except Exception:
            continue
        if pos > 0:
            positions.append(pos)
    return sorted(set(positions))


def _result_is_successful(entry: Dict[str, Any]) -> bool:
    return entry.get("status") == "Success" and entry.get("workspace") and entry.get("design")


def _find_screening_structure_file(workspace_name: str, design_name: str, project_name: str) -> Optional[str]:
    """Locate model_0 PDB for a screening job."""
    if not workspace_name or not design_name or not project_name:
        return None
    filename = f"{workspace_name}_{design_name}_model_0.pdb"
    boltz_results_dir = f"boltz_results_{workspace_name}_{design_name}"
    predictions_dir = f"{workspace_name}_{design_name}"

    candidates = [
        os.path.join(
            "boltzomics_screening_results",
            project_name,
            boltz_results_dir,
            "predictions",
            predictions_dir,
            filename,
        ),
        os.path.join("boltzomics_screening_results", project_name, filename),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def _consistency_label(score: float) -> str:
    if score >= 0.75:
        return "HIGH"
    if score >= 0.5:
        return "MEDIUM"
    return "LOW"


def _best_wt_by_drug(results: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    wt_map: Dict[str, Dict[str, Any]] = {}
    for entry in results:
        if not _result_is_successful(entry):
            continue
        if not _is_wt_name(str(entry.get("protein_name", ""))):
            continue
        drug_name = str(entry.get("drug_name", "") or "")
        current = wt_map.get(drug_name)
        if current is None:
            wt_map[drug_name] = entry
            continue
        current_conf = float(current.get("confidence") or 0.0)
        incoming_conf = float(entry.get("confidence") or 0.0)
        if incoming_conf >= current_conf:
            wt_map[drug_name] = entry
    return wt_map


def _compute_pair_consistency(
    wt_pdb_path: str,
    mutant_pdb_path: str,
    mutation_positions: List[int],
) -> Dict[str, Any]:
    if compare_wt_mutant_quick is None:
        return {"error": "structure_refinement_unavailable"}
    if not mutation_positions:
        return {"error": "no_mutation_positions"}
    if not os.path.exists(wt_pdb_path) or not os.path.exists(mutant_pdb_path):
        return {"error": "missing_pdb"}

    per_pos: List[Dict[str, Any]] = []
    for pos in mutation_positions:
        try:
            comparison = compare_wt_mutant_quick(
                wt_pdb=wt_pdb_path,
                mutant_pdb=mutant_pdb_path,
                mutation_position=int(pos),
                auto_detect_chains=True,
            )
            per_pos.append(comparison)
        except Exception as exc:
            per_pos.append({"error": str(exc), "position": pos})

    valid = [item for item in per_pos if isinstance(item, dict) and "error" not in item]
    if not valid:
        return {"error": "comparison_failed", "per_position": per_pos}

    fingerprint_similarity = sum(float(v.get("fingerprint_similarity", 0.0)) for v in valid) / len(valid)
    disruption = sum(float(v.get("disruption_score", 1.0)) for v in valid) / len(valid)
    consistency_score = max(0.0, min(1.0, 0.7 * fingerprint_similarity + 0.3 * (1.0 - disruption)))

    return {
        "score": consistency_score,
        "label": _consistency_label(consistency_score),
        "fingerprint_similarity": fingerprint_similarity,
        "disruption_score": disruption,
        "num_positions": len(mutation_positions),
        "positions": mutation_positions,
        "per_position": per_pos,
    }


def annotate_results_with_mutation_local_consistency(
    results: List[Dict[str, Any]],
    project_name: str,
    enabled: bool = False,
    cache: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, int]]:
    """
    Annotate result rows with mutation-local consistency metrics.

    Returns:
        (updated_results, cache, summary_counts)
    """
    if cache is None:
        cache = {}

    if not enabled or not results or not project_name:
        return results, cache, {"annotated": 0, "skipped": len(results), "errors": 0}

    updated = [dict(item) for item in results]
    wt_by_drug = _best_wt_by_drug(updated)
    summary = {"annotated": 0, "skipped": 0, "errors": 0}

    for entry in updated:
        if not _result_is_successful(entry):
            summary["skipped"] += 1
            continue
        if _is_wt_name(str(entry.get("protein_name", ""))):
            summary["skipped"] += 1
            continue

        positions = _extract_mutation_positions(str(entry.get("protein_name", "")))
        if not positions:
            summary["skipped"] += 1
            continue

        drug_name = str(entry.get("drug_name", "") or "")
        wt_entry = wt_by_drug.get(drug_name)
        if not wt_entry:
            entry["mutation_local_consistency_reason"] = "missing_wt_reference"
            summary["errors"] += 1
            continue

        wt_ws = str(wt_entry.get("workspace", ""))
        wt_design = str(wt_entry.get("design", ""))
        mut_ws = str(entry.get("workspace", ""))
        mut_design = str(entry.get("design", ""))
        cache_key = (
            f"{project_name}|{wt_ws}|{wt_design}|{mut_ws}|{mut_design}|"
            f"{','.join(str(p) for p in positions)}"
        )

        metrics = cache.get(cache_key)
        if metrics is None:
            wt_pdb = _find_screening_structure_file(wt_ws, wt_design, project_name)
            mut_pdb = _find_screening_structure_file(mut_ws, mut_design, project_name)
            if not wt_pdb or not mut_pdb:
                metrics = {"error": "missing_pdb"}
            else:
                metrics = _compute_pair_consistency(wt_pdb, mut_pdb, positions)
            cache[cache_key] = metrics

        if "error" in metrics:
            entry["mutation_local_consistency_reason"] = str(metrics.get("error"))
            summary["errors"] += 1
            continue

        entry["mutation_local_consistency_score"] = float(metrics["score"])
        entry["mutation_local_consistency_label"] = str(metrics["label"])
        entry["mutation_local_fingerprint_similarity"] = float(metrics["fingerprint_similarity"])
        entry["mutation_local_disruption_score"] = float(metrics["disruption_score"])
        entry["mutation_local_consistency_n_mutations"] = int(metrics["num_positions"])
        entry["mutation_local_consistency_positions"] = ",".join(str(p) for p in metrics["positions"])
        entry["mutation_local_consistency_reason"] = "ok"
        summary["annotated"] += 1

    return updated, cache, summary

