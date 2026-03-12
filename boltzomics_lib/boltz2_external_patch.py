"""
External runtime patch layer for Boltz2 affinity inference.

This module monkeypatches Boltz2 at runtime (without editing the installed package)
to apply a modular affinity post-selection correction and to persist patch metadata
into affinity JSON outputs.
"""

from __future__ import annotations

import functools
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import torch


@dataclass
class ExternalPatchConfig:
    enabled: bool = False
    mode: str = "mutation_aware_v2"
    weight_floor: float = 0.05
    entropy_alpha: float = 0.20
    uncertainty_penalty: float = 0.15
    min_confidence: float = 0.35
    topk: int = 3
    selector_temperature: float = 0.35
    selector_w_iptm: float = 0.35
    selector_w_ligand_iptm: float = 0.25
    selector_w_complex_plddt: float = 0.15
    selector_w_local: float = 0.35
    selector_w_complex_ipde: float = 0.20
    local_contact_cutoff: float = 8.0
    local_pde_weight: float = 0.30
    mutation_positions: List[int] = None
    debug: bool = False


def _env_bool(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def _env_float(name: str, default: float) -> float:
    raw = os.getenv(name)
    if raw is None:
        return float(default)
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float(default)


def _env_int(name: str, default: int) -> int:
    raw = os.getenv(name)
    if raw is None:
        return int(default)
    try:
        return int(raw)
    except (TypeError, ValueError):
        return int(default)


def _env_int_list(name: str) -> List[int]:
    raw = os.getenv(name, "")
    if not raw:
        return []
    out: List[int] = []
    for token in re.split(r"[,\s;:|]+", str(raw).strip()):
        if not token:
            continue
        try:
            value = int(token)
        except (TypeError, ValueError):
            continue
        if value > 0 and value not in out:
            out.append(value)
    return out


def load_config_from_env() -> ExternalPatchConfig:
    return ExternalPatchConfig(
        enabled=_env_bool("BOLTZ2_EXTERNAL_PATCH_ENABLED", False),
        mode=str(os.getenv("BOLTZ2_EXTERNAL_PATCH_MODE", "mutation_aware_v2")).strip().lower(),
        weight_floor=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_WEIGHT_FLOOR", 0.05)),
        entropy_alpha=min(1.0, max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_ENTROPY_ALPHA", 0.20))),
        uncertainty_penalty=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_UNCERTAINTY_PENALTY", 0.15)),
        min_confidence=min(1.0, max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_MIN_CONFIDENCE", 0.35))),
        topk=max(1, _env_int("BOLTZ2_EXTERNAL_PATCH_TOPK", 3)),
        selector_temperature=max(1e-3, _env_float("BOLTZ2_EXTERNAL_PATCH_SELECTOR_TEMPERATURE", 0.35)),
        selector_w_iptm=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_SELECTOR_W_IPTM", 0.35)),
        selector_w_ligand_iptm=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_SELECTOR_W_LIGAND_IPTM", 0.25)),
        selector_w_complex_plddt=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_SELECTOR_W_COMPLEX_PLDDT", 0.15)),
        selector_w_local=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_SELECTOR_W_LOCAL", 0.35)),
        selector_w_complex_ipde=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_SELECTOR_W_COMPLEX_IPDE", 0.20)),
        local_contact_cutoff=max(2.0, _env_float("BOLTZ2_EXTERNAL_PATCH_LOCAL_CONTACT_CUTOFF", 8.0)),
        local_pde_weight=max(0.0, _env_float("BOLTZ2_EXTERNAL_PATCH_LOCAL_PDE_WEIGHT", 0.30)),
        mutation_positions=_env_int_list("BOLTZ2_EXTERNAL_PATCH_MUTATION_POSITIONS"),
        debug=_env_bool("BOLTZ2_EXTERNAL_PATCH_DEBUG", False),
    )


def _to_scalar_float(x: Any, default: float = 0.0) -> float:
    if x is None:
        return float(default)
    if torch.is_tensor(x):
        if x.numel() == 0:
            return float(default)
        return float(torch.nanmean(x.float()).item())
    try:
        return float(x)
    except (TypeError, ValueError):
        return float(default)


def _as_like_tensor(reference: Any, value: float) -> torch.Tensor:
    if torch.is_tensor(reference):
        return torch.full_like(reference, float(value))
    return torch.tensor(float(value), dtype=torch.float32)


def _entropy_binary(prob: float) -> float:
    p = min(1.0 - 1e-8, max(1e-8, float(prob)))
    return -(p * torch.log(torch.tensor(p)) + (1.0 - p) * torch.log(torch.tensor(1.0 - p))).item() / float(torch.log(torch.tensor(2.0)).item())


def _to_1d_float_tensor(
    value: Any,
    *,
    length: Optional[int] = None,
    default: float = 0.0,
    device: Optional[torch.device] = None,
) -> torch.Tensor:
    if torch.is_tensor(value):
        out = value.detach().float().reshape(-1)
        if device is not None:
            out = out.to(device)
    elif isinstance(value, (list, tuple)):
        out = torch.tensor([float(v) for v in value], dtype=torch.float32, device=device)
    else:
        out = torch.tensor([float(_to_scalar_float(value, default))], dtype=torch.float32, device=device)

    if length is None:
        return out
    if out.numel() == length:
        return out
    if out.numel() == 1:
        return out.repeat(length)
    if out.numel() > length:
        return out[:length]
    pad = torch.full((length - out.numel(),), float(default), dtype=out.dtype, device=out.device)
    return torch.cat([out, pad], dim=0)


def _decode_index_vector(value: Any) -> Optional[torch.Tensor]:
    if value is None or (not torch.is_tensor(value)):
        return None
    tensor = value.detach()
    if tensor.dim() == 3:
        tensor = tensor[0]
    if tensor.dim() == 2:
        return torch.argmax(tensor, dim=-1).long()
    if tensor.dim() == 1:
        return tensor.long()
    return None


def _safe_mean(value: torch.Tensor, default: float = 0.0) -> float:
    if value.numel() == 0:
        return float(default)
    return float(value.float().mean().item())


def _build_selector_scores(
    prediction: Dict[str, Any],
    feats: Dict[str, Any],
    sample_atom_coords: torch.Tensor,
    cfg: ExternalPatchConfig,
) -> Dict[str, Any]:
    k = int(sample_atom_coords.shape[0])
    device = sample_atom_coords.device

    iptm = _to_1d_float_tensor(prediction.get("iptm"), length=k, default=0.0, device=device)
    ligand_iptm = _to_1d_float_tensor(
        prediction.get("ligand_iptm"),
        length=k,
        default=float(_safe_mean(iptm, 0.0)),
        device=device,
    )
    complex_plddt = _to_1d_float_tensor(
        prediction.get("complex_plddt"),
        length=k,
        default=0.7,
        device=device,
    )
    complex_ipde = _to_1d_float_tensor(
        prediction.get("complex_ipde"),
        length=k,
        default=1.0,
        device=device,
    )

    token_pad_mask = feats.get("token_pad_mask")
    affinity_token_mask = feats.get("affinity_token_mask")
    mol_type = feats.get("mol_type")
    residue_index = feats.get("residue_index")
    token_to_rep_atom = feats.get("token_to_rep_atom")
    plddt = prediction.get("plddt")
    pde = prediction.get("pde")

    local_scores: List[float] = []
    local_plddt_scores: List[float] = []
    local_ipde_scores: List[float] = []

    if not all(torch.is_tensor(x) for x in [token_pad_mask, affinity_token_mask, mol_type, token_to_rep_atom]):
        local_scores = [float(_safe_mean(complex_plddt, 0.7)) for _ in range(k)]
    else:
        pad = token_pad_mask[0].detach().bool()
        lig_mask = affinity_token_mask[0].detach().bool() & pad
        rec_mask = (mol_type[0].detach() == 0) & pad
        rep_atom_idx = _decode_index_vector(token_to_rep_atom)
        mutation_mask = None
        if torch.is_tensor(residue_index) and cfg.mutation_positions:
            ridx = residue_index[0].detach().long()
            allowed = torch.tensor(cfg.mutation_positions, dtype=ridx.dtype, device=ridx.device)
            mutation_mask = torch.isin(ridx, allowed) & rec_mask

        if rep_atom_idx is None:
            local_scores = [float(_safe_mean(complex_plddt, 0.7)) for _ in range(k)]
        else:
            rep_atom_idx = rep_atom_idx.to(device)
            for i in range(k):
                coords = sample_atom_coords[i]
                tok_coords = coords[rep_atom_idx]
                rec_idx = torch.where(rec_mask)[0]
                lig_idx = torch.where(lig_mask)[0]
                if rec_idx.numel() == 0:
                    local_scores.append(float(complex_plddt[i].item()))
                    local_plddt_scores.append(float(complex_plddt[i].item()))
                    local_ipde_scores.append(float(complex_ipde[i].item()))
                    continue

                target_mask = torch.zeros_like(rec_mask)
                if mutation_mask is not None and int(mutation_mask.sum().item()) > 0:
                    mut_idx = torch.where(mutation_mask)[0]
                    d_mut = torch.cdist(tok_coords[rec_idx], tok_coords[mut_idx])
                    near_mut = d_mut.min(dim=1).values <= float(cfg.local_contact_cutoff)
                    near_mut_idx = rec_idx[near_mut]
                    target_mask[near_mut_idx] = True

                if int(target_mask.sum().item()) == 0:
                    if lig_idx.numel() > 0:
                        d_lig = torch.cdist(tok_coords[rec_idx], tok_coords[lig_idx])
                        near_lig = d_lig.min(dim=1).values <= float(cfg.local_contact_cutoff)
                        near_lig_idx = rec_idx[near_lig]
                        if near_lig_idx.numel() > 0:
                            target_mask[near_lig_idx] = True

                if int(target_mask.sum().item()) == 0:
                    target_mask = rec_mask.clone()

                local_plddt = float(complex_plddt[i].item())
                if torch.is_tensor(plddt):
                    plddt_i = plddt[i] if plddt.dim() >= 2 else plddt
                    if plddt_i.numel() >= target_mask.numel():
                        local_plddt = _safe_mean(plddt_i[: target_mask.numel()][target_mask], default=local_plddt)

                local_ipde = float(complex_ipde[i].item())
                if torch.is_tensor(pde):
                    pde_i = pde[i] if pde.dim() >= 3 else pde
                    if pde_i.dim() == 2 and pde_i.shape[0] >= target_mask.numel() and pde_i.shape[1] >= target_mask.numel():
                        src = target_mask
                        dst = lig_mask if int(lig_mask.sum().item()) > 0 else target_mask
                        local_pairs = pde_i[: target_mask.numel(), : target_mask.numel()][src][:, dst]
                        local_ipde = _safe_mean(local_pairs, default=local_ipde)

                local_pde_penalty = local_ipde / (1.0 + local_ipde)
                local_score = local_plddt - float(cfg.local_pde_weight) * local_pde_penalty
                local_scores.append(float(local_score))
                local_plddt_scores.append(float(local_plddt))
                local_ipde_scores.append(float(local_ipde))

    local_scores_t = _to_1d_float_tensor(local_scores, length=k, default=0.7, device=device)
    complex_ipde_penalty = complex_ipde / (1.0 + complex_ipde)

    selector_scores = (
        float(cfg.selector_w_iptm) * iptm
        + float(cfg.selector_w_ligand_iptm) * ligand_iptm
        + float(cfg.selector_w_complex_plddt) * complex_plddt
        + float(cfg.selector_w_local) * local_scores_t
        - float(cfg.selector_w_complex_ipde) * complex_ipde_penalty
    )
    return {
        "selector_scores": selector_scores,
        "local_scores": local_scores_t,
        "local_plddt_scores": _to_1d_float_tensor(local_plddt_scores, length=k, default=0.0, device=device),
        "local_ipde_scores": _to_1d_float_tensor(local_ipde_scores, length=k, default=0.0, device=device),
        "iptm": iptm,
        "ligand_iptm": ligand_iptm,
        "complex_plddt": complex_plddt,
        "complex_ipde": complex_ipde,
    }


def _run_affinity_at_conformer(
    model: Any,
    feats: Dict[str, Any],
    s_inputs_affinity: torch.Tensor,
    z_affinity: torch.Tensor,
    coords_affinity: torch.Tensor,
) -> Dict[str, float]:
    with torch.autocast("cuda", enabled=False):
        if model.affinity_ensemble:
            out1 = model.affinity_module1(
                s_inputs=s_inputs_affinity.detach(),
                z=z_affinity.detach(),
                x_pred=coords_affinity,
                feats=feats,
                multiplicity=1,
                use_trifast=model.use_trifast,
            )
            out2 = model.affinity_module2(
                s_inputs=s_inputs_affinity.detach(),
                z=z_affinity.detach(),
                x_pred=coords_affinity,
                feats=feats,
                multiplicity=1,
                use_trifast=model.use_trifast,
            )
            v1 = float(out1["affinity_pred_value"].reshape(-1)[0].item())
            v2 = float(out2["affinity_pred_value"].reshape(-1)[0].item())
            p1 = float(torch.sigmoid(out1["affinity_logits_binary"]).reshape(-1)[0].item())
            p2 = float(torch.sigmoid(out2["affinity_logits_binary"]).reshape(-1)[0].item())
            value = (v1 + v2) / 2.0
            prob = (p1 + p2) / 2.0
            if model.affinity_mw_correction:
                model_coef = 1.03525938
                mw_coef = -0.59992683
                bias = 2.83288489
                mw = float((feats["affinity_mw"][0] ** 0.3).reshape(-1)[0].item())
                value = model_coef * value + mw_coef * mw + bias
            return {
                "value": float(value),
                "prob": float(prob),
                "head1_value": float(v1),
                "head1_prob": float(p1),
                "head2_value": float(v2),
                "head2_prob": float(p2),
            }

        out = model.affinity_module(
            s_inputs=s_inputs_affinity.detach(),
            z=z_affinity.detach(),
            x_pred=coords_affinity,
            feats=feats,
            multiplicity=1,
            use_trifast=model.use_trifast,
        )
        value = float(out["affinity_pred_value"].reshape(-1)[0].item())
        prob = float(torch.sigmoid(out["affinity_logits_binary"]).reshape(-1)[0].item())
        return {"value": float(value), "prob": float(prob)}


def _apply_affinity_patch_v2(
    prediction: Dict[str, Any],
    cfg: ExternalPatchConfig,
    model: Any,
    feats: Optional[Dict[str, Any]],
    captured: Optional[Dict[str, Any]],
) -> bool:
    if feats is None or captured is None:
        return False
    if "sample_atom_coords" not in prediction or "iptm" not in prediction:
        return False
    if not torch.is_tensor(prediction.get("sample_atom_coords")):
        return False

    sample_atom_coords = prediction["sample_atom_coords"].detach()
    if sample_atom_coords.dim() != 3:
        return False
    num_samples = int(sample_atom_coords.shape[0])
    if num_samples < 2:
        return False

    z = captured.get("z")
    if not torch.is_tensor(z):
        return False

    pad_token_mask = feats["token_pad_mask"][0]
    rec_mask = (feats["mol_type"][0] == 0) * pad_token_mask
    lig_mask = feats["affinity_token_mask"][0].to(torch.bool) * pad_token_mask
    cross_pair_mask = (
        lig_mask[:, None] * rec_mask[None, :]
        + rec_mask[:, None] * lig_mask[None, :]
        + lig_mask[:, None] * lig_mask[None, :]
    )
    z_affinity = z.detach() * cross_pair_mask[None, :, :, None]

    selector_data = _build_selector_scores(
        prediction=prediction,
        feats=feats,
        sample_atom_coords=sample_atom_coords,
        cfg=cfg,
    )
    selector_scores = selector_data["selector_scores"]
    ranked_idx = torch.argsort(selector_scores, descending=True)
    topk = int(max(1, min(cfg.topk, num_samples)))
    topk_idx = ranked_idx[:topk]
    topk_scores = selector_scores[topk_idx]
    topk_weights = torch.softmax(topk_scores / float(cfg.selector_temperature), dim=0)

    s_inputs_affinity = model.input_embedder(feats, affinity=True)

    affinity_values: List[float] = []
    affinity_probs: List[float] = []
    head1_values: List[float] = []
    head1_probs: List[float] = []
    head2_values: List[float] = []
    head2_probs: List[float] = []
    for idx in topk_idx.tolist():
        coords_affinity = sample_atom_coords[idx][None, None]
        conformer_out = _run_affinity_at_conformer(
            model=model,
            feats=feats,
            s_inputs_affinity=s_inputs_affinity,
            z_affinity=z_affinity,
            coords_affinity=coords_affinity,
        )
        affinity_values.append(float(conformer_out["value"]))
        affinity_probs.append(float(conformer_out["prob"]))
        if "head1_value" in conformer_out:
            head1_values.append(float(conformer_out["head1_value"]))
            head1_probs.append(float(conformer_out["head1_prob"]))
            head2_values.append(float(conformer_out["head2_value"]))
            head2_probs.append(float(conformer_out["head2_prob"]))

    value_tensor = torch.tensor(affinity_values, dtype=torch.float32, device=topk_weights.device)
    prob_tensor = torch.tensor(affinity_probs, dtype=torch.float32, device=topk_weights.device)
    topk_value = float(torch.sum(value_tensor * topk_weights).item())
    topk_prob = float(torch.sum(prob_tensor * topk_weights).item())

    raw_value_t = prediction["affinity_pred_value"]
    raw_prob_t = prediction["affinity_probability_binary"]
    raw_value = _to_scalar_float(raw_value_t, 0.0)
    raw_prob = _to_scalar_float(raw_prob_t, 0.5)
    raw_best_idx = int(torch.argsort(_to_1d_float_tensor(prediction.get("iptm"), length=num_samples), descending=True)[0].item())

    prediction["affinity_pred_value_raw"] = _as_like_tensor(raw_value_t, raw_value)
    prediction["affinity_probability_binary_raw"] = _as_like_tensor(raw_prob_t, raw_prob)
    prediction["affinity_patch_applied"] = _as_like_tensor(raw_value_t, 1.0)
    prediction["affinity_patch_mode"] = cfg.mode
    prediction["affinity_patch_topk_applied"] = _as_like_tensor(raw_value_t, 1.0)
    prediction["affinity_patch_topk_k"] = _as_like_tensor(raw_value_t, float(topk))
    prediction["affinity_patch_best_idx_raw"] = _as_like_tensor(raw_value_t, float(raw_best_idx))
    prediction["affinity_patch_topk_indices"] = topk_idx.detach().cpu()
    prediction["affinity_patch_topk_weights"] = topk_weights.detach().cpu()
    prediction["affinity_patch_selector_scores"] = selector_scores.detach().cpu()
    prediction["affinity_patch_selector_local_scores"] = selector_data["local_scores"].detach().cpu()
    prediction["affinity_patch_selector_local_plddt"] = selector_data["local_plddt_scores"].detach().cpu()
    prediction["affinity_patch_selector_local_ipde"] = selector_data["local_ipde_scores"].detach().cpu()
    prediction["affinity_patch_selector_mutation_positions"] = list(cfg.mutation_positions or [])
    prediction["affinity_patch_topk_pred_values"] = value_tensor.detach().cpu()
    prediction["affinity_patch_topk_probabilities"] = prob_tensor.detach().cpu()
    prediction["affinity_patch_topk_value_std"] = _as_like_tensor(
        raw_value_t,
        float(value_tensor.std(unbiased=False).item()) if value_tensor.numel() > 1 else 0.0,
    )
    prediction["affinity_patch_topk_prob_std"] = _as_like_tensor(
        raw_prob_t,
        float(prob_tensor.std(unbiased=False).item()) if prob_tensor.numel() > 1 else 0.0,
    )
    prediction["affinity_pred_value_patch_base"] = _as_like_tensor(raw_value_t, topk_value)
    prediction["affinity_pred_value"] = _as_like_tensor(raw_value_t, topk_value)
    prediction["affinity_probability_binary"] = _as_like_tensor(raw_prob_t, topk_prob)

    if head1_values and head2_values:
        h1v = torch.tensor(head1_values, dtype=torch.float32, device=topk_weights.device)
        h2v = torch.tensor(head2_values, dtype=torch.float32, device=topk_weights.device)
        h1p = torch.tensor(head1_probs, dtype=torch.float32, device=topk_weights.device)
        h2p = torch.tensor(head2_probs, dtype=torch.float32, device=topk_weights.device)
        prediction["affinity_pred_value1"] = _as_like_tensor(raw_value_t, float(torch.sum(h1v * topk_weights).item()))
        prediction["affinity_probability_binary1"] = _as_like_tensor(raw_prob_t, float(torch.sum(h1p * topk_weights).item()))
        prediction["affinity_pred_value2"] = _as_like_tensor(raw_value_t, float(torch.sum(h2v * topk_weights).item()))
        prediction["affinity_probability_binary2"] = _as_like_tensor(raw_prob_t, float(torch.sum(h2p * topk_weights).item()))

    return True


def _apply_affinity_patch_v1(prediction: Dict[str, Any], cfg: ExternalPatchConfig) -> None:
    if not cfg.enabled:
        return
    if "affinity_pred_value" not in prediction:
        return
    if "affinity_probability_binary" not in prediction:
        return

    raw_value_t = prediction["affinity_pred_value"]
    raw_prob_t = prediction["affinity_probability_binary"]
    raw_value = _to_scalar_float(raw_value_t, 0.0)
    raw_prob = _to_scalar_float(raw_prob_t, 0.5)

    head_values = [raw_value]
    head_probs = [raw_prob]
    for value_key, prob_key in [
        ("affinity_pred_value1", "affinity_probability_binary1"),
        ("affinity_pred_value2", "affinity_probability_binary2"),
    ]:
        if value_key in prediction:
            head_values.append(_to_scalar_float(prediction.get(value_key), raw_value))
            head_probs.append(_to_scalar_float(prediction.get(prob_key), 0.5))

    if len(head_values) == 1:
        prediction["affinity_patch_applied"] = _as_like_tensor(raw_value_t, 1.0)
        prediction["affinity_patch_mode"] = cfg.mode
        prediction["affinity_pred_value_raw"] = _as_like_tensor(raw_value_t, raw_value)
        prediction["affinity_probability_binary_raw"] = _as_like_tensor(raw_prob_t, raw_prob)
        return

    weights = []
    for p in head_probs:
        center_dist = abs(float(p) - 0.5) * 2.0
        entropy = _entropy_binary(p)
        w = max(cfg.weight_floor, center_dist) * (1.0 - cfg.entropy_alpha * entropy)
        weights.append(max(1e-8, float(w)))
    weight_sum = float(sum(weights))
    if weight_sum <= 1e-8:
        weights = [1.0 for _ in weights]
        weight_sum = float(len(weights))

    fused_value = float(sum(v * w for v, w in zip(head_values, weights)) / weight_sum)
    fused_prob = float(sum(p * w for p, w in zip(head_probs, weights)) / weight_sum)

    # Confidence-gated risk adjustment from available Boltz confidence signals.
    iptm = _to_scalar_float(prediction.get("iptm"), raw_prob)
    ligand_iptm = _to_scalar_float(prediction.get("ligand_iptm"), iptm)
    complex_plddt = _to_scalar_float(prediction.get("complex_plddt"), 0.7)
    complex_ipde = _to_scalar_float(prediction.get("complex_ipde"), 0.0)

    gate = (
        torch.sigmoid(torch.tensor((iptm - cfg.min_confidence) * 8.0)).item()
        * torch.sigmoid(torch.tensor((ligand_iptm - cfg.min_confidence) * 8.0)).item()
        * torch.sigmoid(torch.tensor((complex_plddt - 0.70) * 8.0)).item()
        * torch.sigmoid(torch.tensor((0.80 - complex_ipde) * 3.0)).item()
    )
    spread = float(torch.tensor(head_values, dtype=torch.float32).std(unbiased=False).item())
    prob_uncertainty = 1.0 - abs(fused_prob - 0.5) * 2.0
    penalty = cfg.uncertainty_penalty * spread * prob_uncertainty * max(0.0, 1.0 - gate)
    patched_value = fused_value - penalty

    prediction["affinity_pred_value_raw"] = _as_like_tensor(raw_value_t, raw_value)
    prediction["affinity_probability_binary_raw"] = _as_like_tensor(raw_prob_t, raw_prob)
    prediction["affinity_pred_value_patch_base"] = _as_like_tensor(raw_value_t, fused_value)
    prediction["affinity_patch_confidence_gate"] = _as_like_tensor(raw_value_t, gate)
    prediction["affinity_patch_penalty"] = _as_like_tensor(raw_value_t, penalty)
    prediction["affinity_patch_applied"] = _as_like_tensor(raw_value_t, 1.0)
    prediction["affinity_patch_mode"] = cfg.mode
    prediction["affinity_pred_value"] = _as_like_tensor(raw_value_t, patched_value)
    prediction["affinity_probability_binary"] = _as_like_tensor(raw_prob_t, fused_prob)


def _apply_affinity_patch(
    prediction: Dict[str, Any],
    cfg: ExternalPatchConfig,
    *,
    model: Optional[Any] = None,
    feats: Optional[Dict[str, Any]] = None,
    captured: Optional[Dict[str, Any]] = None,
) -> None:
    mode = str(cfg.mode or "").strip().lower()
    if mode in {"mutation_aware_v2", "mutation_local_topk_v1"}:
        if model is not None and _apply_affinity_patch_v2(prediction, cfg, model, feats, captured):
            return
    _apply_affinity_patch_v1(prediction, cfg)


def _to_json_value(value: Any) -> Optional[Any]:
    if torch.is_tensor(value):
        if value.numel() == 0:
            return None
        if value.numel() == 1:
            return float(value.reshape(-1)[0].item())
        return [float(v) for v in value.detach().cpu().reshape(-1).tolist()]
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        return value
    if isinstance(value, (list, tuple)):
        out: List[Any] = []
        for item in value:
            if isinstance(item, (int, float)):
                out.append(float(item))
            else:
                out.append(item)
        return out
    return None


def apply_runtime_patch(cfg: Optional[ExternalPatchConfig] = None) -> bool:
    cfg = cfg or load_config_from_env()
    if not cfg.enabled:
        return False

    import boltz.model.models.boltz2 as boltz2_mod
    from boltz.data.write.writer import BoltzAffinityWriter

    if getattr(boltz2_mod.Boltz2, "_external_patch_applied", False):
        return True

    original_forward = boltz2_mod.Boltz2.forward

    @functools.wraps(original_forward)
    def forward_with_external_patch(self, *args, **kwargs):
        captured: Dict[str, Any] = {}
        hooks = []

        def _capture_pairformer(_module, _inputs, outputs):
            try:
                if isinstance(outputs, tuple) and len(outputs) >= 2:
                    captured["s"] = outputs[0].detach()
                    captured["z"] = outputs[1].detach()
            except Exception:
                pass

        try:
            hooks.append(self.pairformer_module.register_forward_hook(_capture_pairformer))
        except Exception:
            pass
        try:
            if hasattr(self.pairformer_module, "_orig_mod"):
                hooks.append(self.pairformer_module._orig_mod.register_forward_hook(_capture_pairformer))  # noqa: SLF001
        except Exception:
            pass

        try:
            out = original_forward(self, *args, **kwargs)
        finally:
            for hook in hooks:
                try:
                    hook.remove()
                except Exception:
                    pass

        try:
            feats = kwargs.get("feats")
            if feats is None and args:
                if isinstance(args[0], dict):
                    feats = args[0]
            _apply_affinity_patch(
                out,
                cfg,
                model=self,
                feats=feats,
                captured=captured,
            )
        except Exception as exc:
            if cfg.debug:
                print(f"[BOLTZ2-EXTERNAL-PATCH] forward patch failed: {exc}")
        return out

    boltz2_mod.Boltz2.forward = forward_with_external_patch
    boltz2_mod.Boltz2._external_patch_applied = True

    original_write = BoltzAffinityWriter.write_on_batch_end

    @functools.wraps(original_write)
    def write_with_external_patch(self, trainer, pl_module, prediction, batch_indices, batch, batch_idx, dataloader_idx):
        original_write(self, trainer, pl_module, prediction, batch_indices, batch, batch_idx, dataloader_idx)
        try:
            if prediction.get("exception"):
                return
            record_id = batch["record"][0].id
            path = Path(self.output_dir) / record_id / f"affinity_{record_id}.json"
            if not path.exists():
                return
            with path.open("r", encoding="utf-8") as f:
                data = json.load(f)

            extras = {
                "affinity_pred_value_raw": prediction.get("affinity_pred_value_raw"),
                "affinity_probability_binary_raw": prediction.get("affinity_probability_binary_raw"),
                "affinity_pred_value_patch_base": prediction.get("affinity_pred_value_patch_base"),
                "affinity_patch_confidence_gate": prediction.get("affinity_patch_confidence_gate"),
                "affinity_patch_penalty": prediction.get("affinity_patch_penalty"),
                "affinity_patch_applied": prediction.get("affinity_patch_applied"),
                "affinity_patch_mode": prediction.get("affinity_patch_mode"),
                "affinity_patch_topk_applied": prediction.get("affinity_patch_topk_applied"),
                "affinity_patch_topk_k": prediction.get("affinity_patch_topk_k"),
                "affinity_patch_best_idx_raw": prediction.get("affinity_patch_best_idx_raw"),
                "affinity_patch_topk_indices": prediction.get("affinity_patch_topk_indices"),
                "affinity_patch_topk_weights": prediction.get("affinity_patch_topk_weights"),
                "affinity_patch_topk_pred_values": prediction.get("affinity_patch_topk_pred_values"),
                "affinity_patch_topk_probabilities": prediction.get("affinity_patch_topk_probabilities"),
                "affinity_patch_topk_value_std": prediction.get("affinity_patch_topk_value_std"),
                "affinity_patch_topk_prob_std": prediction.get("affinity_patch_topk_prob_std"),
                "affinity_patch_selector_scores": prediction.get("affinity_patch_selector_scores"),
                "affinity_patch_selector_local_scores": prediction.get("affinity_patch_selector_local_scores"),
                "affinity_patch_selector_local_plddt": prediction.get("affinity_patch_selector_local_plddt"),
                "affinity_patch_selector_local_ipde": prediction.get("affinity_patch_selector_local_ipde"),
                "affinity_patch_selector_mutation_positions": prediction.get("affinity_patch_selector_mutation_positions"),
            }
            changed = False
            for key, value in extras.items():
                parsed = _to_json_value(value)
                if parsed is not None:
                    data[key] = parsed
                    changed = True
            if changed:
                with path.open("w", encoding="utf-8") as f:
                    json.dump(data, f, indent=4)
        except Exception as exc:
            if cfg.debug:
                print(f"[BOLTZ2-EXTERNAL-PATCH] writer patch failed: {exc}")

    BoltzAffinityWriter.write_on_batch_end = write_with_external_patch
    BoltzAffinityWriter._external_patch_applied = True
    return True
