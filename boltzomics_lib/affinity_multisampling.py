"""
Affinity multi-sampling helper.

Reuses an existing Boltz structure/pre-affinity output and reruns only the affinity
stage across multiple affinity settings to estimate a range/uncertainty, without
editing the installed Boltz package.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

import utils


@dataclass
class MultiSamplingResult:
    success: bool
    summary_path: Optional[str]
    aggregated_value: Optional[float]
    aggregated_probability: Optional[float]
    selected_setting: Optional[str]
    message: str = ""


def _safe_int(value: Any) -> Optional[int]:
    try:
        if value is None:
            return None
        out = int(str(value).strip())
        return out
    except (TypeError, ValueError):
        return None


def _normalize_refinement_steps(values: Optional[Iterable[Any]]) -> List[int]:
    if not values:
        return []
    out: List[int] = []
    for value in values:
        step = _safe_int(value)
        if step is None or step < 1:
            continue
        if step not in out:
            out.append(step)
    return out


def _profile_label(sampling_steps_affinity: int, diffusion_samples_affinity: int) -> str:
    return f"{int(sampling_steps_affinity)}x{int(diffusion_samples_affinity)}"


def _normalize_profile_list(
    values: Optional[Iterable[Any]],
    diffusion_samples_affinity_base: int,
) -> List[Dict[str, int]]:
    if not values:
        return []
    out: List[Dict[str, int]] = []
    seen: set = set()
    base_diff = max(1, int(diffusion_samples_affinity_base))
    for value in values:
        step: Optional[int] = None
        diff: Optional[int] = None

        if isinstance(value, dict):
            step = _safe_int(
                value.get("sampling_steps_affinity", value.get("sampling_steps", value.get("step")))
            )
            diff = _safe_int(
                value.get("diffusion_samples_affinity", value.get("diffusion_samples", value.get("diff")))
            )
        elif isinstance(value, (tuple, list)):
            if len(value) >= 1:
                step = _safe_int(value[0])
            if len(value) >= 2:
                diff = _safe_int(value[1])
        else:
            token = str(value).strip()
            if token:
                parts = re.split(r"[:/xX]", token, maxsplit=1)
                step = _safe_int(parts[0])
                if len(parts) > 1:
                    diff = _safe_int(parts[1])

        if step is None or step < 1:
            continue
        if diff is None or diff < 1:
            diff = base_diff

        key = (int(step), int(diff))
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "label": _profile_label(step, diff),
                "sampling_steps_affinity": int(step),
                "diffusion_samples_affinity": int(diff),
            }
        )
    return out


def _build_sweep_profiles(
    *,
    profiles: Optional[Iterable[Any]],
    refinement_steps: Optional[Iterable[Any]],
    sampling_steps_affinity_base: int,
    diffusion_samples_affinity_base: int,
) -> List[Dict[str, int]]:
    explicit_profiles = _normalize_profile_list(
        profiles,
        diffusion_samples_affinity_base=int(diffusion_samples_affinity_base),
    )
    if explicit_profiles:
        return explicit_profiles

    # Preferred path: explicit ordered numeric refinement steps from user.
    custom_steps = _normalize_refinement_steps(refinement_steps)
    if custom_steps:
        return [
            {
                "label": _profile_label(step, int(diffusion_samples_affinity_base)),
                "sampling_steps_affinity": int(step),
                "diffusion_samples_affinity": int(diffusion_samples_affinity_base),
            }
            for step in custom_steps
        ]

    # Default: single-profile sweep using current affinity configuration.
    return [
        {
            "label": _profile_label(
                int(sampling_steps_affinity_base),
                int(diffusion_samples_affinity_base),
            ),
            "sampling_steps_affinity": int(sampling_steps_affinity_base),
            "diffusion_samples_affinity": int(diffusion_samples_affinity_base),
        }
    ]


def _read_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(path: Path, data: Dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)


def _safe_float(value: Any) -> Optional[float]:
    try:
        if value is None:
            return None
        v = float(value)
        if np.isfinite(v):
            return v
        return None
    except (TypeError, ValueError):
        return None


def _aggregate(values: np.ndarray, mode: str, weights: Optional[np.ndarray] = None) -> float:
    mode = str(mode or "median").strip().lower()
    if mode == "uncertainty_weighted":
        if weights is None or len(weights) != len(values):
            return float(np.mean(values))
        safe_weights = np.array(weights, dtype=float)
        safe_weights = np.clip(safe_weights, 0.0, None)
        if float(np.sum(safe_weights)) <= 1e-12:
            return float(np.mean(values))
        return float(np.average(values, weights=safe_weights))
    if mode == "mean":
        return float(np.mean(values))
    if mode == "trimmed_mean":
        if len(values) >= 3:
            v = np.sort(values)[1:-1]
            return float(np.mean(v))
        return float(np.mean(values))
    return float(np.median(values))


def _mad_outlier_mask(values: np.ndarray, zmax: float = 3.5) -> np.ndarray:
    if values.size == 0:
        return np.array([], dtype=bool)
    if values.size < 3:
        return np.ones(values.shape[0], dtype=bool)
    median = float(np.median(values))
    abs_dev = np.abs(values - median)
    mad = float(np.median(abs_dev))
    if mad <= 1e-12:
        return np.ones(values.shape[0], dtype=bool)
    robust_z = 0.6745 * abs_dev / mad
    return robust_z <= float(max(1e-6, zmax))


def _bootstrap_ci(
    values: np.ndarray,
    mode: str,
    n_bootstrap: int,
    weights: Optional[np.ndarray] = None,
    alpha: float = 0.05,
) -> Tuple[Optional[float], Optional[float]]:
    if values.size == 0:
        return None, None
    if values.size == 1:
        v = float(values[0])
        return v, v
    n = max(0, int(n_bootstrap))
    if n < 20:
        return None, None
    rng = np.random.default_rng(42)
    estimates = np.empty(n, dtype=float)
    for i in range(n):
        sample_idx = rng.integers(0, values.size, size=values.size)
        sample_values = values[sample_idx]
        sample_weights = None
        if weights is not None and len(weights) == len(values):
            sample_weights = weights[sample_idx]
        estimates[i] = _aggregate(sample_values, mode, weights=sample_weights)
    lo = float(np.quantile(estimates, alpha / 2.0))
    hi = float(np.quantile(estimates, 1.0 - (alpha / 2.0)))
    return lo, hi


def _setting_uncertainty(setting_data: Dict[str, Any]) -> Optional[float]:
    head_values: List[float] = []
    for key, raw_value in (setting_data or {}).items():
        if not re.match(r"^affinity_pred_value\d*$", str(key)):
            continue
        parsed = _safe_float(raw_value)
        if parsed is not None:
            head_values.append(float(parsed))

    if len(head_values) >= 2:
        return float(max(0.0, np.std(np.array(head_values, dtype=float), ddof=1)))

    prob = _safe_float((setting_data or {}).get("affinity_probability_binary"))
    if prob is None:
        return None
    confidence = float(np.clip(abs(prob - 0.5) * 2.0, 0.0, 1.0))
    return float(1.0 - confidence)


def _pick_recommended_setting(
    per_setting: Dict[str, Dict[str, Any]],
    confidence_target: str,
) -> Optional[str]:
    target = str(confidence_target or "balanced").strip().lower()
    ordered = sorted(
        per_setting.items(),
        key=lambda item: (
            int(item[1].get("sampling_steps_affinity", 0) or 0),
            int(item[1].get("diffusion_samples_affinity", 0) or 0),
            str(item[0]),
        ),
    )
    settings = [label for label, _ in ordered]
    if not settings:
        return None
    if target == "high":
        return settings[-1]
    if target == "fast":
        return settings[0]
    return settings[len(settings) // 2]


def run_affinity_multisampling(
    yaml_filepath: str,
    *,
    use_gpu: bool,
    override: bool,
    recycling_steps: int,
    sampling_steps: int,
    diffusion_samples: int,
    max_parallel_samples: int,
    step_scale: float,
    affinity_mw_correction: bool,
    max_msa_seqs: int,
    subsample_msa: bool,
    num_subsampled_msa: int,
    timeout: int,
    use_cached_msa: bool,
    accelerator: str,
    devices: int,
    cuda_visible_devices: Optional[str],
    preprocessing_threads: int,
    use_potentials: bool,
    method: Optional[str],
    external_boltz_patch_enabled: bool,
    external_boltz_patch_mode: str,
    external_boltz_patch_weight_floor: float,
    external_boltz_patch_entropy_alpha: float,
    external_boltz_patch_uncertainty_penalty: float,
    external_boltz_patch_min_confidence: float,
    external_boltz_patch_mutation_positions: Optional[List[int]] = None,
    sampling_steps_affinity_base: int = 300,
    diffusion_samples_affinity_base: int = 7,
    profiles: Optional[Iterable[Any]] = None,
    refinement_steps: Optional[Iterable[Any]] = None,
    aggregate_mode: str = "median",
    apply_aggregate: bool = True,
    confidence_target: str = "balanced",
    early_stop_enabled: bool = True,
    early_stop_min_points: int = 2,
    early_stop_delta: float = 0.02,
    early_stop_std: float = 0.04,
    early_stop_patience: int = 1,
    robust_outlier_filter: bool = True,
    robust_outlier_zmax: float = 3.5,
    bootstrap_samples: int = 300,
) -> MultiSamplingResult:
    yaml_path = Path(yaml_filepath)
    yaml_name = yaml_path.stem
    out_root = yaml_path.parent / f"boltz_results_{yaml_name}" / "predictions" / yaml_name
    affinity_path = out_root / f"affinity_{yaml_name}.json"
    pre_affinity_path = out_root / f"pre_affinity_{yaml_name}.npz"

    def _run_affinity_once(
        *,
        sampling_steps_affinity: int,
        diffusion_samples_affinity: int,
        force_override: bool,
    ) -> Dict[str, Any]:
        utils.run_boltz_prediction(
            yaml_filepath=yaml_filepath,
            use_gpu=use_gpu,
            override=bool(force_override),
            recycling_steps=recycling_steps,
            sampling_steps=sampling_steps,
            diffusion_samples=diffusion_samples,
            max_parallel_samples=max_parallel_samples,
            step_scale=step_scale,
            affinity_mw_correction=affinity_mw_correction,
            max_msa_seqs=max_msa_seqs,
            sampling_steps_affinity=int(sampling_steps_affinity),
            diffusion_samples_affinity=int(diffusion_samples_affinity),
            subsample_msa=subsample_msa,
            num_subsampled_msa=num_subsampled_msa,
            timeout=timeout,
            use_cached_msa=use_cached_msa,
            accelerator=accelerator,
            devices=devices,
            cuda_visible_devices=cuda_visible_devices,
            preprocessing_threads=preprocessing_threads,
            use_potentials=use_potentials,
            method=method,
            external_boltz_patch_enabled=external_boltz_patch_enabled,
            external_boltz_patch_mode=external_boltz_patch_mode,
            external_boltz_patch_weight_floor=external_boltz_patch_weight_floor,
            external_boltz_patch_entropy_alpha=external_boltz_patch_entropy_alpha,
            external_boltz_patch_uncertainty_penalty=external_boltz_patch_uncertainty_penalty,
            external_boltz_patch_min_confidence=external_boltz_patch_min_confidence,
            external_boltz_patch_mutation_positions=external_boltz_patch_mutation_positions,
        )
        if not affinity_path.exists():
            raise FileNotFoundError(f"Missing affinity output after rerun: {affinity_path}")
        return _read_json(affinity_path)

    if not affinity_path.exists():
        return MultiSamplingResult(
            success=False,
            summary_path=None,
            aggregated_value=None,
            aggregated_probability=None,
            selected_setting=None,
            message=f"Missing affinity output: {affinity_path}",
        )

    sweep_profiles = _build_sweep_profiles(
        profiles=profiles,
        refinement_steps=refinement_steps,
        sampling_steps_affinity_base=int(sampling_steps_affinity_base),
        diffusion_samples_affinity_base=int(diffusion_samples_affinity_base),
    )
    requested_settings = [str(profile["label"]) for profile in sweep_profiles]

    if not pre_affinity_path.exists():
        try:
            _run_affinity_once(
                sampling_steps_affinity=int(sampling_steps_affinity_base),
                diffusion_samples_affinity=int(diffusion_samples_affinity_base),
                force_override=True,
            )
        except Exception as exc:
            return MultiSamplingResult(
                success=False,
                summary_path=None,
                aggregated_value=None,
                aggregated_probability=None,
                selected_setting=None,
                message=(
                    "Missing pre-affinity cache required for affinity-only sweep and "
                    f"cache rebuild failed: {exc}"
                ),
            )

    base_data = _read_json(affinity_path)
    initial_steps = int(base_data.get("sampling_steps_affinity", 0) or 0)
    initial_diff = int(base_data.get("diffusion_samples_affinity", 0) or 0)
    current_label = None
    for profile in sweep_profiles:
        label = str(profile["label"])
        if (
            int(profile["sampling_steps_affinity"]) == initial_steps
            and int(profile["diffusion_samples_affinity"]) == initial_diff
        ):
            current_label = label
            break

    per_setting: Dict[str, Dict[str, Any]] = {}
    executed_settings: List[str] = []
    skipped_due_convergence: List[str] = []
    running_trace: List[Dict[str, Any]] = []
    running_values: List[float] = []
    prev_running_agg: Optional[float] = None
    convergence_streak = 0
    early_stop_triggered = False

    for idx, profile in enumerate(sweep_profiles):
        label = str(profile["label"])
        if label == current_label:
            setting_data = dict(base_data)
        else:
            # Force re-run of affinity stage only by removing affinity output;
            # structure outputs remain cached and are reused.
            if affinity_path.exists():
                affinity_path.unlink()
            step_now = int(profile["sampling_steps_affinity"])
            diff_now = int(profile["diffusion_samples_affinity"])
            try:
                setting_data = _run_affinity_once(
                    sampling_steps_affinity=step_now,
                    diffusion_samples_affinity=diff_now,
                    force_override=bool(override),
                )
            except Exception as exc:
                message = str(exc)
                # Recovery path for legacy/incomplete caches missing pre_affinity_*.npz.
                if ("pre_affinity_" in message) or ("FileNotFoundError" in message):
                    setting_data = _run_affinity_once(
                        sampling_steps_affinity=step_now,
                        diffusion_samples_affinity=diff_now,
                        force_override=True,
                    )
                else:
                    raise

        setting_data["sampling_steps_affinity"] = int(profile["sampling_steps_affinity"])
        setting_data["diffusion_samples_affinity"] = int(profile["diffusion_samples_affinity"])
        per_setting[label] = setting_data
        executed_settings.append(label)

        setting_copy = out_root / f"affinity_{yaml_name}_setting_{label}.json"
        _write_json(setting_copy, setting_data)

        # Convergence-aware early stop can shorten long sweeps.
        if early_stop_enabled:
            value_now = _safe_float(setting_data.get("affinity_pred_value"))
            if value_now is not None:
                running_values.append(float(value_now))
            if (
                idx < (len(sweep_profiles) - 1)
                and len(running_values) >= max(2, int(early_stop_min_points))
            ):
                running_arr = np.array(running_values, dtype=float)
                running_agg = _aggregate(running_arr, aggregate_mode)
                running_std = (
                    float(np.std(running_arr, ddof=1))
                    if running_arr.size > 1
                    else 0.0
                )
                delta = (
                    abs(float(running_agg) - float(prev_running_agg))
                    if prev_running_agg is not None
                    else None
                )
                running_trace.append(
                    {
                        "setting": label,
                        "n": int(running_arr.size),
                        "running_aggregate": float(running_agg),
                        "running_std": float(running_std),
                        "delta_vs_prev": float(delta) if delta is not None else None,
                    }
                )
                if (
                    delta is not None
                    and float(delta) <= float(max(0.0, early_stop_delta))
                    and float(running_std) <= float(max(0.0, early_stop_std))
                ):
                    convergence_streak += 1
                else:
                    convergence_streak = 0
                prev_running_agg = float(running_agg)

                if convergence_streak >= max(1, int(early_stop_patience)):
                    early_stop_triggered = True
                    skipped_due_convergence = [
                        str(item["label"]) for item in sweep_profiles[idx + 1 :]
                    ]
                    break

    setting_values_map = {
        label: _safe_float((per_setting.get(label) or {}).get("affinity_pred_value"))
        for label in requested_settings
    }
    setting_prob_map = {
        label: _safe_float((per_setting.get(label) or {}).get("affinity_probability_binary"))
        for label in requested_settings
    }

    valid_labels = [
        label for label in executed_settings
        if setting_values_map.get(label) is not None
    ]
    values: List[float] = [float(setting_values_map[label]) for label in valid_labels]
    probs: List[float] = [
        float(setting_prob_map[label]) if setting_prob_map.get(label) is not None else 0.5
        for label in valid_labels
    ]
    if not values:
        return MultiSamplingResult(
            success=False,
            summary_path=None,
            aggregated_value=None,
            aggregated_probability=None,
            selected_setting=None,
            message="No valid affinity_pred_value across selected settings.",
        )

    values_arr_raw = np.array(values, dtype=float)
    probs_arr_raw = np.array(probs, dtype=float)
    aggregate_labels = list(valid_labels)
    outlier_settings: List[str] = []
    values_arr = values_arr_raw
    probs_arr = probs_arr_raw
    setting_uncertainty_map: Dict[str, Optional[float]] = {
        label: _setting_uncertainty(per_setting.get(label) or {})
        for label in valid_labels
    }
    aggregate_uncertainties: List[Optional[float]] = [setting_uncertainty_map.get(label) for label in aggregate_labels]
    if robust_outlier_filter and values_arr_raw.size >= 3:
        keep_mask = _mad_outlier_mask(values_arr_raw, robust_outlier_zmax)
        kept_n = int(np.sum(keep_mask))
        if 2 <= kept_n < values_arr_raw.size:
            aggregate_labels = [
                label for label, keep in zip(valid_labels, keep_mask) if keep
            ]
            outlier_settings = [
                label for label, keep in zip(valid_labels, keep_mask) if not keep
            ]
            values_arr = values_arr_raw[keep_mask]
            probs_arr = probs_arr_raw[keep_mask]
            aggregate_uncertainties = [
                uncertainty for uncertainty, keep in zip(aggregate_uncertainties, keep_mask) if keep
            ]

    weights_arr: Optional[np.ndarray] = None
    setting_weight_map: Dict[str, Optional[float]] = {label: None for label in requested_settings}
    if str(aggregate_mode or "").strip().lower() == "uncertainty_weighted":
        raw_weights: List[float] = []
        for uncertainty in aggregate_uncertainties:
            if uncertainty is None:
                raw_weights.append(1.0)
            else:
                raw_weights.append(float(1.0 / max(1e-3, float(uncertainty))))
        weights_arr = np.array(raw_weights, dtype=float)
        if float(np.sum(weights_arr)) <= 1e-12:
            weights_arr = np.ones_like(values_arr, dtype=float)
        normalized = weights_arr / float(np.sum(weights_arr))
        for label, weight in zip(aggregate_labels, normalized.tolist()):
            setting_weight_map[str(label)] = float(weight)

    aggregate_value = _aggregate(values_arr, aggregate_mode, weights=weights_arr)
    aggregate_prob = _aggregate(probs_arr, aggregate_mode, weights=weights_arr)
    ci95_low, ci95_high = _bootstrap_ci(
        values_arr,
        aggregate_mode,
        n_bootstrap=bootstrap_samples,
        weights=weights_arr,
        alpha=0.05,
    )

    selected_setting = _pick_recommended_setting(
        {label: per_setting[label] for label in executed_settings if label in per_setting},
        confidence_target=confidence_target,
    )

    summary = {
        "affinity_multisampling_enabled": True,
        "affinity_multisampling_settings": requested_settings,
        "affinity_multisampling_requested_settings": requested_settings,
        "affinity_multisampling_executed_settings": executed_settings,
        "affinity_multisampling_skipped_settings": skipped_due_convergence,
        "affinity_multisampling_requested_n": int(len(requested_settings)),
        "affinity_multisampling_executed_n": int(len(executed_settings)),
        "affinity_multisampling_skipped_n": int(len(skipped_due_convergence)),
        "affinity_multisampling_runtime_saved_fraction_estimate": (
            float(len(skipped_due_convergence) / max(1, len(requested_settings)))
            if len(requested_settings) > 0
            else 0.0
        ),
        "affinity_multisampling_early_stop_enabled": bool(early_stop_enabled),
        "affinity_multisampling_early_stop_triggered": bool(early_stop_triggered),
        "affinity_multisampling_early_stop_min_points": int(max(1, early_stop_min_points)),
        "affinity_multisampling_early_stop_delta": float(max(0.0, early_stop_delta)),
        "affinity_multisampling_early_stop_std": float(max(0.0, early_stop_std)),
        "affinity_multisampling_early_stop_patience": int(max(1, early_stop_patience)),
        "affinity_multisampling_convergence_trace": running_trace,
        "affinity_multisampling_refinement_steps": [
            int(profile["sampling_steps_affinity"]) for profile in sweep_profiles
        ],
        "affinity_multisampling_diffusion_samples": [
            int(profile["diffusion_samples_affinity"]) for profile in sweep_profiles
        ],
        "affinity_multisampling_profiles": [
            {
                "label": str(profile["label"]),
                "sampling_steps_affinity": int(profile["sampling_steps_affinity"]),
                "diffusion_samples_affinity": int(profile["diffusion_samples_affinity"]),
            }
            for profile in sweep_profiles
        ],
        "affinity_multisampling_setting_profiles": {
            str(profile["label"]): {
                "sampling_steps_affinity": int(profile["sampling_steps_affinity"]),
                "diffusion_samples_affinity": int(profile["diffusion_samples_affinity"]),
            }
            for profile in sweep_profiles
        },
        "affinity_multisampling_mode": str(aggregate_mode),
        "affinity_multisampling_apply_aggregate": bool(apply_aggregate),
        "affinity_multisampling_recommended_setting": selected_setting,
        "affinity_multisampling_n": int(len(values_arr_raw)),
        "affinity_multisampling_aggregate_n": int(len(values_arr)),
        "affinity_multisampling_aggregate_settings": aggregate_labels,
        "affinity_multisampling_weighting_mode": (
            "inverse_uncertainty"
            if str(aggregate_mode or "").strip().lower() == "uncertainty_weighted"
            else "none"
        ),
        "affinity_multisampling_setting_uncertainties": {
            str(label): setting_uncertainty_map.get(label)
            for label in requested_settings
        },
        "affinity_multisampling_setting_weights": setting_weight_map,
        "affinity_multisampling_outlier_filter_enabled": bool(robust_outlier_filter),
        "affinity_multisampling_outlier_filter_zmax": float(max(1e-6, robust_outlier_zmax)),
        "affinity_multisampling_outlier_method": "mad_zscore",
        "affinity_multisampling_outlier_settings": outlier_settings,
        "affinity_multisampling_bootstrap_samples": int(max(0, bootstrap_samples)),
        "affinity_multisampling_aggregate": float(aggregate_value),
        "affinity_multisampling_probability_aggregate": float(aggregate_prob),
        "affinity_multisampling_std": float(np.std(values_arr_raw, ddof=1)) if len(values_arr_raw) > 1 else 0.0,
        "affinity_multisampling_min": float(np.min(values_arr_raw)),
        "affinity_multisampling_max": float(np.max(values_arr_raw)),
        "affinity_multisampling_sem": (
            float(np.std(values_arr_raw, ddof=1) / np.sqrt(values_arr_raw.size))
            if values_arr_raw.size > 1
            else 0.0
        ),
        "affinity_multisampling_aggregate_ci95_low": ci95_low,
        "affinity_multisampling_aggregate_ci95_high": ci95_high,
        "affinity_multisampling_aggregate_ci95_width": (
            float(ci95_high - ci95_low)
            if ci95_low is not None and ci95_high is not None
            else None
        ),
        "affinity_multisampling_setting_values": {
            label: setting_values_map.get(label)
            for label in requested_settings
        },
        "affinity_multisampling_setting_probabilities": {
            label: setting_prob_map.get(label)
            for label in requested_settings
        },
    }
    summary_path = out_root / f"affinity_multisampling_{yaml_name}.json"
    _write_json(summary_path, summary)

    final_data = dict(base_data)
    final_data.update(summary)
    final_data["affinity_multisampling_applied"] = bool(apply_aggregate)
    final_data["affinity_pred_value_raw_before_multisampling"] = _safe_float(base_data.get("affinity_pred_value"))
    final_data["affinity_probability_binary_raw_before_multisampling"] = _safe_float(
        base_data.get("affinity_probability_binary")
    )
    if apply_aggregate:
        final_data["affinity_pred_value"] = float(aggregate_value)
        final_data["affinity_probability_binary"] = float(aggregate_prob)
    # Keep primary affinity JSON aligned with multi-sampling metadata.
    _write_json(affinity_path, final_data)

    # Keep per-setting copy of original baseline if available.
    if current_label is not None:
        baseline_copy = out_root / f"affinity_{yaml_name}_setting_{current_label}_baseline.json"
        _write_json(baseline_copy, base_data)

    return MultiSamplingResult(
        success=True,
        summary_path=str(summary_path),
        aggregated_value=float(aggregate_value),
        aggregated_probability=float(aggregate_prob),
        selected_setting=selected_setting,
    )
