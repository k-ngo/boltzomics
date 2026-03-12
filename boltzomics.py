import streamlit as st
import pandas as pd
import numpy as np
import time
import json
import hashlib
import logging
from datetime import datetime
import os
import yaml
from typing import List, Dict, Tuple, Optional, Union, Callable, Any, Set
import plotly.graph_objects as go
from PIL import Image
import plotly.express as px
import shutil
import re
import math
import glob
from tempfile import NamedTemporaryFile
import copy
from collections import OrderedDict
import uuid
import subprocess
import sys
from contextlib import contextmanager

logger = logging.getLogger(__name__)

MODULES_DIR = os.path.join(os.path.dirname(__file__), "boltzomics_lib")
if MODULES_DIR not in sys.path:
    sys.path.insert(0, MODULES_DIR)


try:
    import utils
    from utils import parse_protein_chains
except ImportError:
    utils = None
    parse_protein_chains = None

try:
    from screening_job_manager import ScreeningJobManager, ScreeningJob
except ImportError:
    ScreeningJobManager = None
    ScreeningJob = None

try:
    from msa_cache import MSACache, get_msa_cache, should_use_msa_cache
except ImportError:
    MSACache = None
    get_msa_cache = None
    should_use_msa_cache = None

try:
    from pharmacogenomic_analysis import (
        PharmacogenomicAnalyzer,
        MutationEffect,
        ConfidenceLevel,
        analyze_screening_results,
        quick_delta_analysis
    )
except ImportError:
    PharmacogenomicAnalyzer = None
    MutationEffect = None
    ConfidenceLevel = None
    analyze_screening_results = None
    quick_delta_analysis = None

try:
    from structure_refinement import (
        quick_interface_check,
        compare_wt_mutant_quick,
        analyze_binding_interface,
        generate_interaction_fingerprint,
        refine_and_analyze,
        check_openmm_available
    )
except ImportError:
    quick_interface_check = None
    compare_wt_mutant_quick = None
    analyze_binding_interface = None
    generate_interaction_fingerprint = None
    refine_and_analyze = None
    check_openmm_available = None

try:
    from mutation_analysis import (
        analyze_mutation_impact,
        analyze_mutation_properties,
        analyze_binding_proximity,
        quick_proximity_check,
        quick_property_check
    )
except ImportError:
    analyze_mutation_impact = None
    analyze_mutation_properties = None
    analyze_binding_proximity = None
    quick_proximity_check = None
    quick_property_check = None

try:
    from mutation_local_consistency import annotate_results_with_mutation_local_consistency
except ImportError:
    annotate_results_with_mutation_local_consistency = None

try:
    from affinity_multisampling import run_affinity_multisampling
except ImportError:
    run_affinity_multisampling = None

try:
    import styles
except ImportError:
    styles = None
from drug_screening_input import (
    parse_mutations,
    apply_mutations_to_sequence,
    generate_mutant_name_from_text,
    verify_mutations_with_wt_residues,
    parse_fasta_sequences,
    parse_smiles_list,
    validate_smiles,
    validate_ccd_code,
    validate_protein_sequence,
    display_mutation_discovery_section,
    display_binding_pocket_section,
    display_ptm_section,
    display_protein_drug_filter_section,
    should_evaluate_protein_drug_pair,
    calculate_filtered_job_count
)

# Import project management functions
from project_management import (
    get_project_list,
    load_project_data,
    delete_project,
    save_screening_results,
    rename_results_in_project
)

# Import visualization functions
from drug_screening_visualization import create_visualizations, display_structure_only_3d_viewer, deduplicate_results

# Constants
ESTIMATED_TIME_PER_JOB = 300  # 5 minutes per job in seconds
RESULTS_DIR = "boltzomics_screening_results"
JOB_STATE_DIR = os.path.join(RESULTS_DIR, "_job_state")
USE_SCREENING_JOB_QUEUE = ScreeningJobManager is not None

JOB_MANAGER_SESSION_STATE_KEY = "_screening_job_manager"
JOB_MANAGER_EXECUTOR_STATE_KEY = "_screening_job_manager_executor_registered"
JOB_MANAGER_INIT_ERROR_KEY = "_screening_job_manager_init_error"
QUEUE_STATUS_REFRESH_INTERVAL_SECONDS = 1.0  # Auto-refresh job queue UI every 1s while jobs are active
QUEUE_NEXT_REFRESH_STATE_KEY = "_screening_queue_next_refresh"


def get_job_manager() -> Optional[ScreeningJobManager]:
    """Return a cached job manager instance that survives Streamlit reruns."""
    if not USE_SCREENING_JOB_QUEUE or ScreeningJobManager is None:
        return None

    try:
        session_state = st.session_state
    except Exception:
        session_state = None
    if session_state is None:
        # Fallback for non-Streamlit contexts (e.g., tests)
        try:
            return ScreeningJobManager(JOB_STATE_DIR)
        except Exception:
            return None

    manager = session_state.get(JOB_MANAGER_SESSION_STATE_KEY)
    if manager is None:
        try:
            manager = ScreeningJobManager(JOB_STATE_DIR)
        except Exception as exc:
            session_state[JOB_MANAGER_INIT_ERROR_KEY] = str(exc)
            session_state[JOB_MANAGER_SESSION_STATE_KEY] = None
            return None
        session_state[JOB_MANAGER_SESSION_STATE_KEY] = manager
    return manager


def _format_duration(seconds: Optional[float]) -> Optional[str]:
    if seconds is None:
        return None
    if seconds < 0:
        seconds = 0
    total_seconds = int(round(seconds))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    parts = []
    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    if secs or not parts:
        parts.append(f"{secs}s")
    return " ".join(parts)


def _normalize_for_hash(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _normalize_for_hash(value[k]) for k in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_normalize_for_hash(v) for v in value]
    if isinstance(value, set):
        normalized = [_normalize_for_hash(v) for v in value]
        return sorted(normalized, key=lambda item: json.dumps(item, sort_keys=True))
    return value


def _is_wt_name(name: str) -> bool:
    upper = str(name).strip().upper()
    return (
        upper in {"WT", "WILDTYPE", "WILD_TYPE"}
        or upper.endswith("_WT")
        or upper.endswith("-WT")
    )


def _extract_mutation_positions_from_label(label: str) -> List[int]:
    """
    Parse residue indices from mutation labels (e.g., "K262R", "Y652A-F656A").
    """
    if not label:
        return []
    if _is_wt_name(label):
        return []
    matches = re.findall(r"[A-Za-z\*](\d+)[A-Za-z\*]", str(label))
    out: List[int] = []
    for token in matches:
        try:
            value = int(token)
        except (TypeError, ValueError):
            continue
        if value > 0 and value not in out:
            out.append(value)
    return out


def _order_protein_sequences_for_screening(
    protein_sequences: List[Tuple[str, str]]
) -> List[Tuple[str, str]]:
    """
    Keep input order but move WT-like entries first.
    """
    indexed = list(enumerate(protein_sequences))
    indexed.sort(key=lambda item: (0 if _is_wt_name(item[1][0]) else 1, item[0]))
    return [entry for _, entry in indexed]


def _parse_positive_int_list_from_text(raw_text: str) -> List[int]:
    if not raw_text:
        return []
    tokens = re.split(r"[\s,;]+", str(raw_text).strip())
    out: List[int] = []
    for token in tokens:
        t = token.strip()
        if not t:
            continue
        try:
            value = int(t)
        except ValueError:
            continue
        if value < 1:
            continue
        if value not in out:
            out.append(value)
    return out


def _parse_affinity_multisampling_profiles_from_text(
    raw_text: str,
    default_diffusion_samples: int,
) -> List[Dict[str, int]]:
    if not raw_text:
        return []
    tokens = re.split(r"[\s,;]+", str(raw_text).strip())
    out: List[Dict[str, int]] = []
    seen: set = set()
    base_diff = max(1, int(default_diffusion_samples))
    for token in tokens:
        t = str(token).strip()
        if not t:
            continue
        parts = re.split(r"[:/xX]", t, maxsplit=1)
        try:
            step = int(parts[0].strip())
        except (TypeError, ValueError):
            continue
        if step < 1:
            continue
        diff = base_diff
        if len(parts) > 1 and parts[1].strip():
            try:
                diff = int(parts[1].strip())
            except (TypeError, ValueError):
                continue
        if diff < 1:
            continue
        key = (int(step), int(diff))
        if key in seen:
            continue
        seen.add(key)
        out.append(
            {
                "sampling_steps_affinity": int(step),
                "diffusion_samples_affinity": int(diff),
            }
        )
    return out


def _auto_diffusion_samples_for_affinity_step(step: int) -> int:
    # Anchor to historical presets: 200->5, 300->7, 400->9.
    return max(1, int(round(5.0 + (float(step) - 200.0) / 50.0)))


def _prediction_reproducibility_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Select only prediction-relevant parameters for stable job/output identifiers.
    """
    relevant_keys = [
        "recycling_steps",
        "sampling_steps",
        "diffusion_samples",
        "max_parallel_samples",
        "step_scale",
        "affinity_mw_correction",
        "affinity_consensus_enabled",
        "affinity_consensus_mode",
        "affinity_consensus_weight_floor",
        "affinity_consensus_entropy_alpha",
        "affinity_multisampling_enabled",
        "affinity_multisampling_profiles",
        "affinity_multisampling_settings",
        "affinity_multisampling_refinement_steps",
        "affinity_multisampling_aggregate_mode",
        "affinity_multisampling_apply_aggregate",
        "affinity_multisampling_early_stop_enabled",
        "affinity_multisampling_early_stop_min_points",
        "affinity_multisampling_early_stop_delta",
        "affinity_multisampling_early_stop_std",
        "affinity_multisampling_early_stop_patience",
        "affinity_multisampling_robust_outlier_filter",
        "affinity_multisampling_robust_outlier_zmax",
        "affinity_multisampling_bootstrap_samples",
        "external_boltz_patch_enabled",
        "external_boltz_patch_mode",
        "external_boltz_patch_weight_floor",
        "external_boltz_patch_entropy_alpha",
        "external_boltz_patch_uncertainty_penalty",
        "external_boltz_patch_min_confidence",
        "max_msa_seqs",
        "sampling_steps_affinity",
        "diffusion_samples_affinity",
        "subsample_msa",
        "num_subsampled_msa",
        "use_potentials",
        "method",
        "mutation_steering_config",
        "template_cif_path",
        "binding_pocket_constraints",
        "cofactor_info",
        "ptm_modifications",
        "structure_only",
    ]
    selected: Dict[str, Any] = {}
    for key in relevant_keys:
        if key in params:
            selected[key] = params.get(key)
    return selected


def _derive_stable_job_identifiers(
    project_name: str,
    protein_name: str,
    protein_sequence: str,
    drug_name: str,
    smiles: str,
    structure_only: bool,
    params: Dict[str, Any],
) -> Tuple[str, str, str]:
    """
    Build deterministic identifiers for workspace, design, and queue job id.
    """
    design_name = _sanitize_design_name(protein_name, None if structure_only else drug_name)
    signature_payload = {
        "project_name": project_name,
        "protein_name": protein_name,
        "protein_sequence": protein_sequence,
        "drug_name": "" if structure_only else drug_name,
        "smiles": "" if structure_only else smiles,
        "design_name": design_name,
        "params": _prediction_reproducibility_params(params),
    }
    normalized = _normalize_for_hash(signature_payload)
    signature_raw = json.dumps(normalized, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    digest = hashlib.sha1(signature_raw.encode("utf-8")).hexdigest()
    workspace_name = f"screening_{digest[:12]}"
    job_id = f"{project_name}:{digest[:20]}"
    return workspace_name, design_name, job_id


def _find_existing_results_by_yaml_name(
    project_name: str,
    yaml_name: str,
    structure_only: bool = False,
) -> Optional[str]:
    """
    Fast path for deterministic output IDs before directory-wide scanning.
    """
    if not yaml_name:
        return None
    project_dir = os.path.join(RESULTS_DIR, project_name)
    yaml_filepath = os.path.join(project_dir, f"{yaml_name}.yaml")
    if not os.path.exists(yaml_filepath):
        return None
    is_valid, _ = validate_boltz_results(yaml_filepath, structure_only=structure_only)
    if not is_valid:
        return None
    return yaml_name


def discover_gpu_devices() -> List[Tuple[str, str]]:
    """
    Discover available CUDA GPUs for user-facing selection.

    Returns:
        List of (gpu_index, gpu_name) tuples.
    """
    devices: List[Tuple[str, str]] = []

    # Prefer torch for reliable CUDA device names if available.
    try:
        import torch

        if torch.cuda.is_available():
            count = torch.cuda.device_count()
            for idx in range(count):
                devices.append((str(idx), torch.cuda.get_device_name(idx)))
            if devices:
                return devices
    except Exception:
        pass

    # Fallback to nvidia-smi when torch is unavailable/incomplete.
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        if result.returncode == 0 and result.stdout.strip():
            for line in result.stdout.strip().splitlines():
                parts = [p.strip() for p in line.split(",", 1)]
                if len(parts) == 2:
                    devices.append((parts[0], parts[1]))
    except Exception:
        pass

    return devices

# =============================================================================
# MSA CACHING INTEGRATION
# =============================================================================
# The MSA (Multiple Sequence Alignment) caching system dramatically accelerates
# pharmacogenomic screening by reusing MSA computations across related predictions.
#
# Scientific Justification:
# - MSA captures evolutionary covariance between residue positions
# - Point mutations don't change protein family membership
# - Same homologous sequences are found for WT and mutants
# - Boltz-2 uses MSA for co-evolutionary signal (family-level information)
#
# Performance Impact:
# - MSA generation: ~30-60 seconds per query
# - For N mutants × M drugs: reduces N×M MSA generations to just 1
# - Typical time savings: 95-99% for mutation screening workflows
# =============================================================================

def get_or_create_msa_cache_for_screening(
    wt_sequence: str,
    project_name: str,
    first_drug_smiles: str,
    boltz_params: Dict[str, Any],
    enable_msa_cache: bool = True
) -> Tuple[Optional[str], bool]:
    """
    Get cached MSA for a protein sequence, or prepare for MSA generation.

    This function implements the "generate once, reuse many" strategy for MSA
    caching in mutation screening workflows.

    Args:
        wt_sequence: Wild-type protein sequence
        project_name: Name of the screening project
        first_drug_smiles: SMILES of first drug (used if we need to generate MSA)
        boltz_params: Boltz prediction parameters
        enable_msa_cache: Whether to use MSA caching

    Returns:
        Tuple of (msa_path or None, should_use_cached_msa)
    """
    if not enable_msa_cache or get_msa_cache is None:
        return None, False

    cache = get_msa_cache()

    # Check if MSA already exists for this protein
    existing_msa = cache.get_cached_msa(wt_sequence)
    if existing_msa:
        return existing_msa, True

    # MSA not cached - return None so caller generates it
    # The MSA will be cached after the first prediction completes
    return None, False


def cache_msa_after_prediction(
    protein_sequence: str,
    yaml_filepath: str,
    enable_msa_cache: bool = True
) -> Optional[str]:
    """
    Cache the MSA generated by a Boltz prediction for future reuse.

    Call this after a successful Boltz prediction to cache the MSA.

    Args:
        protein_sequence: The protein sequence that was predicted
        yaml_filepath: Path to the YAML file used for prediction
        enable_msa_cache: Whether MSA caching is enabled

    Returns:
        Path to cached MSA, or None if caching failed/disabled
    """
    if not enable_msa_cache or get_msa_cache is None:
        return None

    try:
        cache = get_msa_cache()

        # Construct path to Boltz output directory
        yaml_dir = os.path.dirname(yaml_filepath)
        yaml_name = os.path.splitext(os.path.basename(yaml_filepath))[0]
        boltz_output_dir = os.path.join(yaml_dir, f"boltz_results_{yaml_name}")

        if os.path.exists(boltz_output_dir):
            return cache.cache_msa_from_boltz_output(
                protein_sequence,
                boltz_output_dir,
                source_info={"yaml_file": yaml_filepath}
            )
    except Exception as e:
        # MSA caching is an optimization - don't fail the prediction if caching fails
        logger.warning("Failed to cache MSA: %s", e, exc_info=True)

    return None


def get_msa_for_mutant(
    wt_sequence: str,
    mutant_sequence: str,
    mutations: List[Tuple[str, int, str, str]] = None,
    enable_msa_cache: bool = True
) -> Tuple[Optional[str], bool]:
    """
    Get MSA for a mutant sequence, leveraging WT MSA cache.

    For point mutations, we can reuse the WT MSA since the evolutionary
    context is essentially identical.

    Args:
        wt_sequence: Wild-type protein sequence
        mutant_sequence: Mutant protein sequence
        mutations: List of mutations as (chain_id, position, wt_res, mut_res)
        enable_msa_cache: Whether MSA caching is enabled

    Returns:
        Tuple of (msa_path or None, should_use_cached_msa)
    """
    if not enable_msa_cache or get_msa_cache is None:
        return None, False

    cache = get_msa_cache()

    # Prefer a mutant-specific cache entry if present.
    mutant_msa = cache.get_cached_msa(mutant_sequence)
    if mutant_msa:
        return mutant_msa, True

    # Fall back to WT cache; when available, attempt to create a mutant-specific
    # query row for best sequence consistency while reusing homolog rows.
    wt_msa = cache.get_cached_msa(wt_sequence)
    if wt_msa:
        try:
            mutant_specific_msa = cache.create_mutant_msa(
                wt_sequence=wt_sequence,
                mutant_sequence=mutant_sequence,
                mutations=mutations or [],
            )
            if mutant_specific_msa:
                return mutant_specific_msa, True
        except Exception as e:
            # MSA caching is an optimization path; safe fallback to WT cache.
            logger.warning("Failed to create mutant-specific MSA, falling back to WT MSA: %s", e, exc_info=True)
        return wt_msa, True

    return None, False


def resolve_msa_cache_for_prediction(
    protein_sequence: str,
    wt_sequence: Optional[str] = None,
    enable_msa_cache: bool = True,
) -> Tuple[Optional[str], bool]:
    """
    Resolve whether a cached MSA can be used for this prediction.

    Returns:
        Tuple of (msa_path, use_cached_msa)
    """
    if not enable_msa_cache:
        return None, False

    try:
        if wt_sequence and wt_sequence != protein_sequence:
            msa_path, use_cached = get_msa_for_mutant(
                wt_sequence=wt_sequence,
                mutant_sequence=protein_sequence,
                enable_msa_cache=enable_msa_cache,
            )
            if msa_path and use_cached:
                return msa_path, True

        if should_use_msa_cache is not None:
            use_cached, msa_path = should_use_msa_cache(
                protein_sequence=protein_sequence,
                wt_sequence=None,
                enable_cache=enable_msa_cache,
            )
            if use_cached and msa_path:
                return msa_path, True
    except Exception as e:
        logger.warning("Failed to resolve cached MSA, will generate fresh: %s", e, exc_info=True)

    return None, False


def _extract_mutation_positions_from_name(protein_name: str) -> List[int]:
    """Extract residue numbers from mutation labels like Y652A-F656A."""
    if not protein_name:
        return []
    matches = re.findall(r"[A-Za-z](\d+)[A-Za-z]", str(protein_name))
    positions: List[int] = []
    for m in matches:
        try:
            pos = int(m)
        except Exception:
            continue
        if pos > 0:
            positions.append(pos)
    return sorted(set(positions))


@contextmanager
def _boltz_job_file_lock(
    lock_path: str,
    timeout_seconds: float = 900.0,
    stale_seconds: float = 21600.0,
):
    """
    Simple cross-process lock using atomic lock-file creation.
    Prevents two workers/sessions from writing the same boltz_results_<yaml> path.
    """
    start = time.time()
    while True:
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(f"pid={os.getpid()} started={time.time():.3f}\n")
            break
        except FileExistsError:
            try:
                mtime = os.path.getmtime(lock_path)
                if (time.time() - mtime) > float(stale_seconds):
                    os.remove(lock_path)
                    continue
            except FileNotFoundError:
                continue
            if (time.time() - start) >= float(timeout_seconds):
                raise TimeoutError(f"Timed out waiting for job lock: {lock_path}")
            time.sleep(0.5)
    try:
        yield
    finally:
        try:
            os.remove(lock_path)
        except FileNotFoundError:
            pass


def _resolve_mutation_steering_constraints(
    protein_name: str,
    protein_sequence: str,
    base_constraints: Optional[Dict[str, Any]],
    structure_only: bool,
    mutation_steering_config: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """
    Build mutation-neighborhood pocket constraints, preserving manual constraints when present.
    """
    if structure_only:
        return base_constraints
    if not mutation_steering_config or not mutation_steering_config.get("enabled", False):
        return base_constraints
    if not mutation_steering_config.get("mutation_mode_active", False):
        return base_constraints
    if _is_wt_name(protein_name) and not mutation_steering_config.get("apply_to_wt", False):
        return base_constraints

    mutation_positions = _extract_mutation_positions_from_name(protein_name)
    if not mutation_positions:
        return base_constraints

    is_valid, _, chains_dict, _ = validate_protein_sequence(protein_sequence)
    if not is_valid or not chains_dict:
        return base_constraints

    chain_starts = mutation_steering_config.get("chain_starts") or {}
    neighborhood = int(mutation_steering_config.get("window", 3))
    neighborhood = max(0, neighborhood)
    steering_max_distance = float(mutation_steering_config.get("max_distance", 6.0))

    contacts_set: Set[Tuple[str, int]] = set()
    if base_constraints and base_constraints.get("contacts"):
        for chain_id, residue_idx in base_constraints.get("contacts", []):
            try:
                contacts_set.add((str(chain_id), int(residue_idx)))
            except Exception:
                continue

    for residue_number in mutation_positions:
        for chain_id, chain_seq in chains_dict.items():
            chain_start = int(chain_starts.get(chain_id, 1))
            chain_end = chain_start + len(chain_seq) - 1
            if not (chain_start <= residue_number <= chain_end):
                continue
            local_residue = residue_number - chain_start + 1
            for idx in range(local_residue - neighborhood, local_residue + neighborhood + 1):
                if 1 <= idx <= len(chain_seq):
                    contacts_set.add((chain_id, idx))

    if not contacts_set:
        return base_constraints

    sorted_contacts = sorted(contacts_set, key=lambda item: (item[0], item[1]))
    binder = (base_constraints or {}).get("binder", "X")
    max_distance = float((base_constraints or {}).get("max_distance", steering_max_distance))
    mode = (base_constraints or {}).get("mode", "mutation_steering")
    if base_constraints and base_constraints.get("contacts"):
        mode = f"{mode}+mutation_steering"

    return {
        "binder": binder,
        "contacts": [[chain_id, residue_idx] for chain_id, residue_idx in sorted_contacts],
        "max_distance": max_distance,
        "mode": mode,
        "source": "mutation_neighborhood",
        "mutation_positions": mutation_positions,
        "neighborhood_window": neighborhood,
    }


def load_css():
    """Loads custom CSS for styling the application."""
    if styles:
        for css in [styles.header_css, styles.bubble_css, styles.button_css,
                   styles.dialog_and_plot_css, styles.modify_streamlit_style]:
            st.markdown(css, unsafe_allow_html=True)
    st.markdown("""
    <style>
    [data-testid='stSidebar'] [data-testid='stTooltipHoverTarget'] svg {
        fill: white !important;
        color: white !important;
    }
    </style>
    """, unsafe_allow_html=True)

def create_screening_boltz_yaml(
    workspace_name,
    design_name,
    protein_sequence,
    ligand_smiles,
    project_name,
    binding_pocket_constraints=None,
    cofactor_info=None,
    template_cif_path=None,
    structure_only=False,
    ptm_modifications=None,
    msa_path=None,
    target_directory: Optional[str] = None,
):
    """Create a YAML file for Boltz prediction in the project-specific directory.

    Args:
        workspace_name: Name of the workspace/project
        design_name: Name of this specific design/prediction
        protein_sequence: Protein amino acid sequence (chains separated by ':')
        ligand_smiles: SMILES string for the ligand
        project_name: Name of the screening project
        binding_pocket_constraints: Optional binding pocket constraints
        cofactor_info: Optional cofactor information
        template_cif_path: Optional path to template CIF file
        structure_only: If True, skip affinity prediction
        ptm_modifications: Optional post-translational modifications
        msa_path: Optional path to cached MSA file. When provided, Boltz will use
                  this pre-computed MSA instead of querying the MSA server, which
                  dramatically speeds up mutation screening (95%+ time savings).

    Returns:
        Path to the created YAML file
    """
    # Create project directory
    project_dir = target_directory or os.path.join("boltzomics_screening_results", project_name)
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)

    # Create filename
    filename = f"{workspace_name}_{design_name}.yaml"
    filepath = os.path.join(project_dir, filename)

    # Parse protein sequence into chains (with optional MSA path for caching)
    try:
        protein_chains = parse_protein_chains(protein_sequence, msa_path=msa_path)
    except ValueError as e:
        raise ValueError(f"Invalid protein sequence format: {str(e)}")
    
    # Add PTM modifications to protein chains if provided
    if ptm_modifications and ptm_modifications.get('modifications'):
        modifications = ptm_modifications.get('modifications', [])
        for mod in modifications:
            chain_id = mod.get('chain_id', 'A')
            position = mod.get('position')
            ccd = mod.get('ccd')
            
            if chain_id and position and ccd:
                # Find the corresponding protein chain and add modification
                for chain in protein_chains:
                    if chain.get('protein', {}).get('id') == chain_id:
                        if 'modifications' not in chain['protein']:
                            chain['protein']['modifications'] = []
                        chain['protein']['modifications'].append({
                            'position': position,
                            'ccd': ccd
                        })
                        break
    
    yaml_content = {"sequences": protein_chains}
    if not structure_only:
        yaml_content["sequences"] += [
            {"ligand": {"id": "X", "smiles": ligand_smiles}}
        ]
    
    # Add co-factors if provided
    if cofactor_info and isinstance(cofactor_info, list) and len(cofactor_info) > 0:
        for i, cofactor in enumerate(cofactor_info):
            if i >= 4:  # Limit to 4 cofactors
                break
            if cofactor and (cofactor.get('smiles') or cofactor.get('ccd')):
                # Use chain IDs T, U, V, W
                chain_id = chr(ord('T') + i)
                cofactor_entry = {
                    "ligand": {
                        "id": chain_id
                    }
                }
                
                # Add either SMILES or CCD code
                if cofactor.get('smiles'):
                    cofactor_entry["ligand"]["smiles"] = cofactor['smiles']
                elif cofactor.get('ccd'):
                    cofactor_entry["ligand"]["ccd"] = cofactor['ccd']
                
                # Add co-factor to sequences
                yaml_content["sequences"].append(cofactor_entry)

    # Add templates section if template_cif_path is provided
    if template_cif_path:
        yaml_content["templates"] = [{"cif": os.path.abspath(template_cif_path)}]
    
    # Add constraints if provided and valid
    if binding_pocket_constraints and binding_pocket_constraints.get('contacts'):
        contacts = []
        for c in binding_pocket_constraints.get('contacts', []):
            if len(c) >= 2:
                # Try to cast residue index to int if possible
                try:
                    res_idx = int(c[1])
                except (ValueError, TypeError):
                    res_idx = c[1]
                contacts.append([c[0], res_idx])
        pocket_constraint = {
            "pocket": {
                "binder": binding_pocket_constraints.get('binder', 'X'),
                "contacts": contacts,
                "max_distance": float(binding_pocket_constraints.get('max_distance', 5.0))
            }
        }
        yaml_content_copy = copy.deepcopy(yaml_content)
        yaml_content_copy["constraints"] = [pocket_constraint]
        constraints = yaml_content_copy.pop("constraints")
        with open(filepath, 'w') as f:
            yaml.dump(yaml_content_copy, f, default_flow_style=False)
            f.write("constraints:\n")
            for constraint in constraints:
                f.write("  - pocket:\n")
                f.write(f"      binder: {constraint['pocket']['binder']}\n")
                contacts_str = yaml.dump(constraint['pocket']['contacts'], default_flow_style=True).strip()
                f.write(f"      contacts: {contacts_str}\n")
                f.write(f"      max_distance: {constraint['pocket']['max_distance']}\n")
            if not structure_only:
                f.write("properties:\n")
                f.write("  - affinity:\n")
                f.write("      binder: X\n")
    else:
        with open(filepath, 'w') as f:
            yaml.dump(yaml_content, f, default_flow_style=False)
            if not structure_only:
                f.write("properties:\n")
                f.write("  - affinity:\n")
                f.write("      binder: X\n")
    
    return filepath

def validate_boltz_results(yaml_filepath, structure_only=False):
    """
    Validate that Boltz results are complete and valid.
    Args:
        yaml_filepath: Path to the YAML file used for prediction
        structure_only: If True, only check for confidence file
    Returns:
        tuple: (is_valid, error_message)
    """
    try:
        # Construct the path to the results files
        yaml_dir = os.path.dirname(yaml_filepath)
        yaml_filename = os.path.basename(yaml_filepath)
        yaml_name = os.path.splitext(yaml_filename)[0]
        # Path to the confidence results JSON file
        confidence_results_path = os.path.join(yaml_dir, f"boltz_results_{yaml_name}", "predictions", yaml_name, f"confidence_{yaml_name}_model_0.json")
        if not structure_only:
            # Path to the affinity results JSON file
            affinity_results_path = os.path.join(yaml_dir, f"boltz_results_{yaml_name}", "predictions", yaml_name, f"affinity_{yaml_name}.json")
            # Check if both files exist
            if not os.path.exists(affinity_results_path):
                return False, f"Affinity results file not found: {affinity_results_path}"
            if not os.path.exists(confidence_results_path):
                return False, f"Confidence results file not found: {confidence_results_path}"
            # Check if files are not empty
            if os.path.getsize(affinity_results_path) == 0:
                return False, f"Affinity results file is empty: {affinity_results_path}"
            if os.path.getsize(confidence_results_path) == 0:
                return False, f"Confidence results file is empty: {confidence_results_path}"
            # Try to parse the files to ensure they contain valid JSON
            try:
                with open(affinity_results_path, 'r') as f:
                    affinity_data = json.load(f)
                # Check for required fields in affinity results
                required_affinity_fields = ["affinity_pred_value", "affinity_probability_binary"]
                for field in required_affinity_fields:
                    if field not in affinity_data:
                        return False, f"Missing required field '{field}' in affinity results"
                    if affinity_data[field] is None:
                        return False, f"Required field '{field}' is None in affinity results"
                    # Check if the values are numeric
                    if not isinstance(affinity_data[field], (int, float)):
                        return False, f"Required field '{field}' is not numeric in affinity results"
            except json.JSONDecodeError as e:
                return False, f"Invalid JSON in affinity results file: {str(e)}"
            except Exception as e:
                return False, f"Error reading affinity results file: {str(e)}"
        else:
            # Structure-only: only check confidence file
            if not os.path.exists(confidence_results_path):
                return False, f"Confidence results file not found: {confidence_results_path}"
            if os.path.getsize(confidence_results_path) == 0:
                return False, f"Confidence results file is empty: {confidence_results_path}"
        # Try to parse the confidence file to ensure it contains valid JSON
        try:
            with open(confidence_results_path, 'r') as f:
                confidence_data = json.load(f)
            # Check for required fields in confidence results
            required_confidence_fields = ["confidence_score", "ptm", "iptm", "complex_plddt"]
            for field in required_confidence_fields:
                if field not in confidence_data:
                    return False, f"Missing required field '{field}' in confidence results"
                if confidence_data[field] is None:
                    return False, f"Required field '{field}' is None in confidence results"
                # Check if the values are numeric
                if not isinstance(confidence_data[field], (int, float)):
                    return False, f"Required field '{field}' is not numeric in confidence results"
        except json.JSONDecodeError as e:
            return False, f"Invalid JSON in confidence results file: {str(e)}"
        except Exception as e:
            return False, f"Error reading confidence results file: {str(e)}"
        return True, "Results validation successful"
    except Exception as e:
        return False, f"Error during results validation: {str(e)}"

def run_boltz_with_retry(
    workspace_name,
    design_name,
    protein_sequence,
    ligand_smiles,
    project_name,
    protein_display_name,
    ligand_display_name,
    use_gpu=True,
    binding_pocket_constraints=None,
    override=False,
    recycling_steps=3,
    sampling_steps=200,
    diffusion_samples=1,
    max_parallel_samples=5,
    step_scale=1.638,
    affinity_mw_correction=False,
    external_boltz_patch_enabled: bool = False,
    external_boltz_patch_mode: str = "mutation_aware_v2",
    external_boltz_patch_weight_floor: float = 0.05,
    external_boltz_patch_entropy_alpha: float = 0.20,
    external_boltz_patch_uncertainty_penalty: float = 0.15,
    external_boltz_patch_min_confidence: float = 0.35,
    affinity_multisampling_enabled: bool = False,
    affinity_multisampling_profiles: Optional[List[Dict[str, int]]] = None,
    affinity_multisampling_settings: Optional[List[str]] = None,
    affinity_multisampling_refinement_steps: Optional[List[int]] = None,
    affinity_multisampling_aggregate_mode: str = "median",
    affinity_multisampling_apply_aggregate: bool = True,
    affinity_multisampling_early_stop_enabled: bool = True,
    affinity_multisampling_early_stop_min_points: int = 2,
    affinity_multisampling_early_stop_delta: float = 0.02,
    affinity_multisampling_early_stop_std: float = 0.04,
    affinity_multisampling_early_stop_patience: int = 1,
    affinity_multisampling_robust_outlier_filter: bool = True,
    affinity_multisampling_robust_outlier_zmax: float = 3.5,
    affinity_multisampling_bootstrap_samples: int = 300,
    confidence_target: str = "balanced",
    max_msa_seqs=8192,
    sampling_steps_affinity=200,
    diffusion_samples_affinity=5,
    cofactor_info=None,
    enable_retries=True,
    max_retry_attempts=2,
    retry_delay_base=5,
    subsample_msa=False,
    num_subsampled_msa=1024,
    template_cif_path=None,
    structure_only=False,
    ptm_modifications=None,
    prediction_timeout_seconds=300,
    accelerator: str = "gpu",
    devices: int = 1,
    cuda_visible_devices: Optional[str] = None,
    preprocessing_threads: int = 1,
    use_potentials: bool = False,
    method: Optional[str] = None,
    emit_streamlit_feedback=True,
    status_callback: Optional[Callable[[str, Dict[str, Union[str, int, float]]], None]] = None,
    msa_path: Optional[str] = None,
    use_cached_msa: bool = False,
    enable_msa_cache: bool = True,
    execution_directory: Optional[str] = None,
):
    """Run Boltz workflow with retry logic and validation.

    Args:
        ... (existing params)
        msa_path: Optional path to cached MSA file. When provided, uses pre-computed
                  MSA instead of querying MSA server, dramatically speeding up
                  mutation screening.
        use_cached_msa: If True and msa_path is provided, skip --use_msa_server flag.
        enable_msa_cache: If True, cache newly generated MSA files after successful runs.
    """
    last_error = None
    total_attempts = max_retry_attempts + 1 if enable_retries else 1
    current_msa_path = msa_path
    current_use_cached_msa = bool(use_cached_msa)
    cache_fallback_triggered = False
    force_override_rebuild = False
    force_clean_rebuild = False
    oom_fallback_applied = False
    current_max_parallel_samples = int(max_parallel_samples)
    current_devices = int(devices)
    patch_mutation_positions = (
        _extract_mutation_positions_from_label(protein_display_name)
        if external_boltz_patch_enabled
        else []
    )

    def notify(event: str, payload: Dict[str, Union[str, int, float]]) -> None:
        if status_callback:
            try:
                status_callback(event, payload)
            except Exception:
                pass

    for attempt in range(total_attempts):
        try:
            yaml_filepath = create_screening_boltz_yaml(
                workspace_name,
                design_name,
                protein_sequence,
                ligand_smiles,
                project_name,
                binding_pocket_constraints,
                cofactor_info,
                template_cif_path,
                structure_only,
                ptm_modifications,
                msa_path=current_msa_path,  # MSA caching support
                target_directory=execution_directory,
            )
            yaml_dir = os.path.dirname(yaml_filepath)
            yaml_name = os.path.splitext(os.path.basename(yaml_filepath))[0]
            boltz_output_dir = os.path.join(yaml_dir, f"boltz_results_{yaml_name}")
            lock_path = os.path.join(yaml_dir, f".{yaml_name}.boltz.lock")
            pre_affinity_path = os.path.join(
                yaml_dir,
                f"boltz_results_{yaml_name}",
                "predictions",
                yaml_name,
                f"pre_affinity_{yaml_name}.npz",
            )
            if force_clean_rebuild and os.path.exists(boltz_output_dir):
                try:
                    shutil.rmtree(boltz_output_dir)
                except Exception:
                    pass
                force_clean_rebuild = False
            if (
                (not structure_only)
                and (not override)
                and (not force_override_rebuild)
                and os.path.exists(os.path.dirname(pre_affinity_path))
                and (not os.path.exists(pre_affinity_path))
            ):
                force_override_rebuild = True
                notify(
                    "pre_affinity_rebuild",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                    },
                )
            with _boltz_job_file_lock(lock_path):
                utils.run_boltz_prediction(
                    yaml_filepath=yaml_filepath,
                    use_gpu=use_gpu,
                    override=(override or force_override_rebuild),
                    recycling_steps=recycling_steps,
                    sampling_steps=sampling_steps,
                    diffusion_samples=diffusion_samples,
                    max_parallel_samples=current_max_parallel_samples,
                    step_scale=step_scale,
                    affinity_mw_correction=affinity_mw_correction,
                    external_boltz_patch_enabled=external_boltz_patch_enabled,
                    external_boltz_patch_mode=external_boltz_patch_mode,
                    external_boltz_patch_weight_floor=external_boltz_patch_weight_floor,
                    external_boltz_patch_entropy_alpha=external_boltz_patch_entropy_alpha,
                    external_boltz_patch_uncertainty_penalty=external_boltz_patch_uncertainty_penalty,
                    external_boltz_patch_min_confidence=external_boltz_patch_min_confidence,
                    external_boltz_patch_mutation_positions=patch_mutation_positions,
                    max_msa_seqs=max_msa_seqs,
                    sampling_steps_affinity=sampling_steps_affinity,
                    diffusion_samples_affinity=diffusion_samples_affinity,
                    subsample_msa=subsample_msa,
                    num_subsampled_msa=num_subsampled_msa,
                    timeout=prediction_timeout_seconds,
                    use_cached_msa=current_use_cached_msa,  # Skip MSA server when using cache
                    accelerator=accelerator,
                    devices=current_devices,
                    cuda_visible_devices=cuda_visible_devices,
                    preprocessing_threads=preprocessing_threads,
                    use_potentials=use_potentials,
                    method=method,
                )
            is_valid, validation_error = validate_boltz_results(yaml_filepath, structure_only=structure_only)
            if not is_valid:
                raise Exception(f"Results validation failed: {validation_error}")
            results = utils.parse_boltz_results(yaml_filepath, structure_only=structure_only)
            if results is None:
                raise Exception("Failed to parse Boltz results")

            # Multi-sampling reruns affinity-only stages that depend on pre_affinity cache.
            # Even in a fresh project folder, an earlier failed attempt for the same job can
            # leave a partial output directory missing this artifact.
            if affinity_multisampling_enabled and (not structure_only):
                yaml_dir = os.path.dirname(yaml_filepath)
                yaml_name = os.path.splitext(os.path.basename(yaml_filepath))[0]
                pre_affinity_path = os.path.join(
                    yaml_dir,
                    f"boltz_results_{yaml_name}",
                    "predictions",
                    yaml_name,
                    f"pre_affinity_{yaml_name}.npz",
                )
                if not os.path.exists(pre_affinity_path):
                    notify(
                        "pre_affinity_rebuild_before_multisampling",
                        {
                            "attempt": attempt + 1,
                            "protein": protein_display_name,
                            "ligand": ligand_display_name,
                        },
                    )
                    with _boltz_job_file_lock(lock_path):
                        utils.run_boltz_prediction(
                            yaml_filepath=yaml_filepath,
                            use_gpu=use_gpu,
                            override=True,
                            recycling_steps=recycling_steps,
                            sampling_steps=sampling_steps,
                            diffusion_samples=diffusion_samples,
                            max_parallel_samples=current_max_parallel_samples,
                            step_scale=step_scale,
                            affinity_mw_correction=affinity_mw_correction,
                            external_boltz_patch_enabled=external_boltz_patch_enabled,
                            external_boltz_patch_mode=external_boltz_patch_mode,
                            external_boltz_patch_weight_floor=external_boltz_patch_weight_floor,
                            external_boltz_patch_entropy_alpha=external_boltz_patch_entropy_alpha,
                            external_boltz_patch_uncertainty_penalty=external_boltz_patch_uncertainty_penalty,
                            external_boltz_patch_min_confidence=external_boltz_patch_min_confidence,
                            external_boltz_patch_mutation_positions=patch_mutation_positions,
                            max_msa_seqs=max_msa_seqs,
                            sampling_steps_affinity=sampling_steps_affinity,
                            diffusion_samples_affinity=diffusion_samples_affinity,
                            subsample_msa=subsample_msa,
                            num_subsampled_msa=num_subsampled_msa,
                            timeout=prediction_timeout_seconds,
                            use_cached_msa=current_use_cached_msa,
                            accelerator=accelerator,
                            devices=current_devices,
                            cuda_visible_devices=cuda_visible_devices,
                            preprocessing_threads=preprocessing_threads,
                            use_potentials=use_potentials,
                            method=method,
                        )
                    is_valid, validation_error = validate_boltz_results(yaml_filepath, structure_only=structure_only)
                    if not is_valid:
                        raise Exception(f"Results validation failed after pre-affinity rebuild: {validation_error}")
                    results = utils.parse_boltz_results(yaml_filepath, structure_only=structure_only)
                    if results is None:
                        raise Exception("Failed to parse Boltz results after pre-affinity rebuild")

            if (
                affinity_multisampling_enabled
                and (not structure_only)
                and run_affinity_multisampling is not None
            ):
                with _boltz_job_file_lock(lock_path):
                    multi_result = run_affinity_multisampling(
                        yaml_filepath=yaml_filepath,
                        use_gpu=use_gpu,
                        override=False,  # preserve prior behavior: reuse existing structure outputs
                        recycling_steps=recycling_steps,
                        sampling_steps=sampling_steps,
                        diffusion_samples=diffusion_samples,
                        max_parallel_samples=current_max_parallel_samples,
                        step_scale=step_scale,
                        affinity_mw_correction=affinity_mw_correction,
                        max_msa_seqs=max_msa_seqs,
                        subsample_msa=subsample_msa,
                        num_subsampled_msa=num_subsampled_msa,
                        timeout=prediction_timeout_seconds,
                        use_cached_msa=current_use_cached_msa,
                        accelerator=accelerator,
                        devices=current_devices,
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
                        external_boltz_patch_mutation_positions=patch_mutation_positions,
                        sampling_steps_affinity_base=sampling_steps_affinity,
                        diffusion_samples_affinity_base=diffusion_samples_affinity,
                        profiles=affinity_multisampling_profiles,
                        refinement_steps=affinity_multisampling_refinement_steps,
                        aggregate_mode=affinity_multisampling_aggregate_mode,
                        apply_aggregate=affinity_multisampling_apply_aggregate,
                        confidence_target=confidence_target,
                        early_stop_enabled=affinity_multisampling_early_stop_enabled,
                        early_stop_min_points=affinity_multisampling_early_stop_min_points,
                        early_stop_delta=affinity_multisampling_early_stop_delta,
                        early_stop_std=affinity_multisampling_early_stop_std,
                        early_stop_patience=affinity_multisampling_early_stop_patience,
                        robust_outlier_filter=affinity_multisampling_robust_outlier_filter,
                        robust_outlier_zmax=affinity_multisampling_robust_outlier_zmax,
                        bootstrap_samples=affinity_multisampling_bootstrap_samples,
                    )
                if multi_result.success:
                    results = utils.parse_boltz_results(yaml_filepath, structure_only=structure_only) or results

            # Cache newly generated MSA so subsequent mutant/drug jobs can reuse it.
            if enable_msa_cache and not current_use_cached_msa:
                cache_msa_after_prediction(
                    protein_sequence=protein_sequence,
                    yaml_filepath=yaml_filepath,
                    enable_msa_cache=enable_msa_cache,
                )

            # Track MSA source in results for post-hoc analysis
            if results and isinstance(results, dict):
                if cache_fallback_triggered:
                    results["msa_source"] = "fresh_after_fallback"
                    results["msa_cache_fallback"] = True
                    results["msa_cache_fallback_error"] = (last_error or "")[:500]
                elif current_use_cached_msa:
                    results["msa_source"] = "cached"
                else:
                    results["msa_source"] = "fresh"

            if cache_fallback_triggered and emit_streamlit_feedback:
                st.info(
                    f"Note: {protein_display_name} + {ligand_display_name} succeeded after "
                    f"MSA cache fallback (original cached-MSA error: {(last_error or '')[:100]})"
                )

            if attempt > 0:
                notify(
                    "success",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                    },
                )
                if emit_streamlit_feedback:
                    st.success(f"Job succeeded for {protein_display_name} + {ligand_display_name} after {attempt+1} attempt(s)!")
            return results, True, None
        except Exception as e:
            last_error = str(e)
            last_error_lower = last_error.lower()
            oom_failure = (
                ("out of memory" in last_error_lower)
                or ("ran out of memory" in last_error_lower)
                or ("number of failed examples" in last_error_lower)
            )
            missing_pre_affinity_cache = (
                "pre_affinity_" in last_error_lower
                and "not found" in last_error_lower
            )
            if (
                oom_failure
                and attempt < max_retry_attempts
                and not oom_fallback_applied
            ):
                oom_fallback_applied = True
                force_override_rebuild = True
                force_clean_rebuild = True
                current_max_parallel_samples = 1
                if str(accelerator).lower() == "gpu":
                    current_devices = 1
                notify(
                    "oom_fallback_retry",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                        "error": last_error[:200],
                    },
                )
                if emit_streamlit_feedback:
                    st.warning(
                        "Boltz reported GPU memory pressure. Retrying with safer memory settings "
                        "(max_parallel_samples=1, devices=1)."
                    )
                continue
            if (
                (not structure_only)
                and missing_pre_affinity_cache
                and (not override)
                and (not force_override_rebuild)
                and attempt < max_retry_attempts
            ):
                force_override_rebuild = True
                force_clean_rebuild = True
                notify(
                    "pre_affinity_rebuild_retry",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                        "error": last_error[:200],
                    },
                )
                if emit_streamlit_feedback:
                    st.warning(
                        "Detected incomplete cached outputs; retrying this job from a clean output folder."
                    )
                continue

            # Robust fallback: some cached-MSA runs can complete structure/confidence
            # yet miss affinity output; retry with fresh MSA generation.
            missing_affinity_output = (
                "Affinity results file not found" in last_error
                or "affinity_" in last_error.lower() and "not found" in last_error.lower()
            )
            cached_input_processing_failure = "input processing failure" in last_error.lower()
            if (
                current_use_cached_msa
                and (missing_affinity_output or cached_input_processing_failure)
                and attempt < max_retry_attempts
            ):
                current_use_cached_msa = False
                current_msa_path = None
                cache_fallback_triggered = True
                notify(
                    "msa_cache_fallback",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                        "error": last_error[:200],
                    },
                )
                if emit_streamlit_feedback:
                    st.warning(
                        "Cached MSA run produced incomplete affinity output; retrying with fresh MSA generation."
                    )
                continue

            if enable_retries and attempt < max_retry_attempts:
                delay = retry_delay_base * (2 ** attempt)
                notify(
                    "retry",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                        "delay": delay,
                        "error": str(e)[:200],
                    },
                )
                if emit_streamlit_feedback:
                    st.warning(f"Attempt {attempt + 1} failed for {protein_display_name} + {ligand_display_name}: {str(e)[:100]}. Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                notify(
                    "failure",
                    {
                        "attempt": attempt + 1,
                        "protein": protein_display_name,
                        "ligand": ligand_display_name,
                        "error": str(e)[:200],
                        "attempts": total_attempts,
                    },
                )
                if emit_streamlit_feedback:
                    if enable_retries:
                        st.error(f"All {total_attempts} attempts failed for {protein_display_name} + {ligand_display_name}: {str(e)[:100]}")
                    else:
                        st.error(f"Prediction failed for {protein_display_name} + {ligand_display_name}: {str(e)[:100]}")
    if cache_fallback_triggered:
        notify(
            "msa_cache_fallback_exhausted",
            {
                "protein": protein_display_name,
                "ligand": ligand_display_name,
                "error": (last_error or "")[:200],
            },
        )
    return None, False, last_error


def _sanitize_design_name(*parts: str) -> str:
    raw = "_".join([p for p in parts if p])
    raw = raw.strip()
    if not raw:
        raw = "screening_job"
    sanitized = re.sub(r"[^a-zA-Z0-9_\-]+", "_", raw)
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    if not sanitized:
        sanitized = f"screening_{uuid.uuid4().hex[:6]}"
    return sanitized[:120]


def _generate_workspace_name(prefix: str = "screening") -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{prefix}_{timestamp}_{uuid.uuid4().hex[:6]}"


def _build_boltz_command(yaml_filename: str, params: Dict[str, Any]) -> str:
    external_patch_enabled = bool(params.get("external_boltz_patch_enabled", False))
    if external_patch_enabled:
        patch_cli = os.path.join(MODULES_DIR, "boltz2_patched_cli.py")
        cmd: List[str] = ["python", patch_cli, "predict", yaml_filename]
    else:
        cmd = ["boltz", "predict", yaml_filename]
    if not params.get("use_cached_msa", False):
        cmd.append("--use_msa_server")
    cmd.extend(["--output_format", "pdb"])
    accelerator = str(params.get("accelerator", "gpu")).lower()
    cmd.extend(["--accelerator", accelerator])
    cmd.extend(["--devices", str(int(params.get("devices", 1)))])
    cmd.extend(["--preprocessing-threads", str(int(params.get("preprocessing_threads", 1)))])
    if params.get("override"):
        cmd.append("--override")
    cmd.extend(["--recycling_steps", str(int(params.get("recycling_steps", 3)))])
    cmd.extend(["--sampling_steps", str(int(params.get("sampling_steps", 200)))])
    cmd.extend(["--diffusion_samples", str(int(params.get("diffusion_samples", 1)))])
    cmd.extend(["--max_parallel_samples", str(int(params.get("max_parallel_samples", 5)))])
    cmd.extend(["--step_scale", str(float(params.get("step_scale", 1.638)))])
    if params.get("affinity_mw_correction"):
        cmd.append("--affinity_mw_correction")
    if params.get("use_potentials"):
        cmd.append("--use_potentials")
    method = params.get("method")
    if method:
        cmd.extend(["--method", str(method)])
    cmd.extend(["--max_msa_seqs", str(int(params.get("max_msa_seqs", 8192)))])
    cmd.extend(["--sampling_steps_affinity", str(int(params.get("sampling_steps_affinity", 200)))])
    cmd.extend(["--diffusion_samples_affinity", str(int(params.get("diffusion_samples_affinity", 5)))])
    if params.get("subsample_msa"):
        cmd.append("--subsample_msa")
        cmd.extend(["--num_subsampled_msa", str(int(params.get("num_subsampled_msa", 1024)))])
    cmd_str = " ".join(cmd)
    if external_patch_enabled:
        env_parts = [
            "BOLTZ2_EXTERNAL_PATCH_ENABLED=1",
            f"BOLTZ2_EXTERNAL_PATCH_MODE={params.get('external_boltz_patch_mode', 'mutation_aware_v2')}",
            f"BOLTZ2_EXTERNAL_PATCH_WEIGHT_FLOOR={float(params.get('external_boltz_patch_weight_floor', 0.05))}",
            f"BOLTZ2_EXTERNAL_PATCH_ENTROPY_ALPHA={float(params.get('external_boltz_patch_entropy_alpha', 0.20))}",
            f"BOLTZ2_EXTERNAL_PATCH_UNCERTAINTY_PENALTY={float(params.get('external_boltz_patch_uncertainty_penalty', 0.15))}",
            f"BOLTZ2_EXTERNAL_PATCH_MIN_CONFIDENCE={float(params.get('external_boltz_patch_min_confidence', 0.35))}",
        ]
        cmd_str = f"{' '.join(env_parts)} {cmd_str}"
    cuda_visible_devices = params.get("cuda_visible_devices")
    if accelerator == "gpu" and cuda_visible_devices not in (None, "", "auto"):
        cmd_str = f"CUDA_VISIBLE_DEVICES={cuda_visible_devices} {cmd_str}"
    return cmd_str


def _build_boltz_batch_command(input_path: str, params: Dict[str, Any], use_msa_server: bool) -> str:
    external_patch_enabled = bool(params.get("external_boltz_patch_enabled", False))
    if external_patch_enabled:
        patch_cli = os.path.join(MODULES_DIR, "boltz2_patched_cli.py")
        cmd: List[str] = ["python", patch_cli, "predict", input_path]
    else:
        cmd = ["boltz", "predict", input_path]
    if use_msa_server:
        cmd.append("--use_msa_server")
    cmd.extend(["--output_format", "pdb"])
    accelerator = str(params.get("accelerator", "gpu")).lower()
    cmd.extend(["--accelerator", accelerator])
    cmd.extend(["--devices", str(int(params.get("devices", 1)))])
    cmd.extend(["--preprocessing-threads", str(int(params.get("preprocessing_threads", 1)))])
    if params.get("override"):
        cmd.append("--override")
    cmd.extend(["--recycling_steps", str(int(params.get("recycling_steps", 3)))])
    cmd.extend(["--sampling_steps", str(int(params.get("sampling_steps", 200)))])
    cmd.extend(["--diffusion_samples", str(int(params.get("diffusion_samples", 1)))])
    cmd.extend(["--max_parallel_samples", str(int(params.get("max_parallel_samples", 5)))])
    cmd.extend(["--step_scale", str(float(params.get("step_scale", 1.638)))])
    if params.get("affinity_mw_correction"):
        cmd.append("--affinity_mw_correction")
    if params.get("use_potentials"):
        cmd.append("--use_potentials")
    method = params.get("method")
    if method:
        cmd.extend(["--method", str(method)])
    cmd.extend(["--max_msa_seqs", str(int(params.get("max_msa_seqs", 8192)))])
    cmd.extend(["--sampling_steps_affinity", str(int(params.get("sampling_steps_affinity", 200)))])
    cmd.extend(["--diffusion_samples_affinity", str(int(params.get("diffusion_samples_affinity", 5)))])
    if params.get("subsample_msa"):
        cmd.append("--subsample_msa")
        cmd.extend(["--num_subsampled_msa", str(int(params.get("num_subsampled_msa", 1024)))])
    cmd_str = " ".join(cmd)
    if external_patch_enabled:
        env_parts = [
            "BOLTZ2_EXTERNAL_PATCH_ENABLED=1",
            f"BOLTZ2_EXTERNAL_PATCH_MODE={params.get('external_boltz_patch_mode', 'mutation_aware_v2')}",
            f"BOLTZ2_EXTERNAL_PATCH_WEIGHT_FLOOR={float(params.get('external_boltz_patch_weight_floor', 0.05))}",
            f"BOLTZ2_EXTERNAL_PATCH_ENTROPY_ALPHA={float(params.get('external_boltz_patch_entropy_alpha', 0.20))}",
            f"BOLTZ2_EXTERNAL_PATCH_UNCERTAINTY_PENALTY={float(params.get('external_boltz_patch_uncertainty_penalty', 0.15))}",
            f"BOLTZ2_EXTERNAL_PATCH_MIN_CONFIDENCE={float(params.get('external_boltz_patch_min_confidence', 0.35))}",
        ]
        cmd_str = f"{' '.join(env_parts)} {cmd_str}"
    cuda_visible_devices = params.get("cuda_visible_devices")
    if accelerator == "gpu" and cuda_visible_devices not in (None, "", "auto"):
        cmd_str = f"CUDA_VISIBLE_DEVICES={cuda_visible_devices} {cmd_str}"
    return cmd_str


def _extract_workspace_design(yaml_name: str) -> Tuple[str, str]:
    if not yaml_name:
        return "screening", "screening"
    return yaml_name, yaml_name


def _to_float_or_none(value: Any) -> Optional[float]:
    try:
        if value is None:
            return None
        v = float(value)
        if np.isfinite(v):
            return v
        return None
    except (TypeError, ValueError):
        return None


def _compute_affinity_consensus(
    boltz_results: Dict[str, Any],
    params: Dict[str, Any],
) -> Dict[str, Any]:
    # Multi-sampling stores both raw single-run and aggregated values.
    # Keep raw as the selected table/CSV value; keep aggregate for plotting/analysis.
    if bool(boltz_results.get("affinity_multisampling_applied", False)):
        aggregated_value = _to_float_or_none(
            boltz_results.get("affinity_multisampling_aggregate", boltz_results.get("affinity_pred_value"))
        )
        aggregated_prob = _to_float_or_none(
            boltz_results.get(
                "affinity_multisampling_probability_aggregate",
                boltz_results.get("affinity_probability_binary"),
            )
        )
        selected_value = _to_float_or_none(
            boltz_results.get(
                "affinity_pred_value_raw_before_multisampling",
                boltz_results.get("affinity_pred_value"),
            )
        )
        selected_prob = _to_float_or_none(
            boltz_results.get(
                "affinity_probability_binary_raw_before_multisampling",
                boltz_results.get("affinity_probability_binary"),
            )
        )
        if aggregated_value is None:
            aggregated_value = selected_value
        if aggregated_prob is None:
            aggregated_prob = selected_prob
        if selected_value is None:
            selected_value = 0.0
        if selected_prob is None:
            selected_prob = 0.0
        return {
            "selected_affinity_pred_value": float(selected_value),
            "selected_affinity_probability": float(selected_prob),
            "raw_affinity_pred_value": _to_float_or_none(boltz_results.get("affinity_pred_value_raw_before_multisampling")),
            "raw_affinity_probability": _to_float_or_none(boltz_results.get("affinity_probability_binary_raw_before_multisampling")),
            "consensus_enabled": True,
            "consensus_mode": "multisampling_aggregate",
            "consensus_n_heads": int(boltz_results.get("affinity_multisampling_n", 1) or 1),
            "consensus_pred_value": float(aggregated_value if aggregated_value is not None else selected_value),
            "consensus_probability": float(aggregated_prob if aggregated_prob is not None else selected_prob),
            "consensus_std": _to_float_or_none(boltz_results.get("affinity_multisampling_std")),
            "consensus_head_values": {},
            "consensus_head_probabilities": {},
        }

    raw_value = _to_float_or_none(boltz_results.get("affinity_pred_value"))
    raw_prob = _to_float_or_none(boltz_results.get("affinity_probability_binary"))

    heads = [
        ("main", raw_value, raw_prob),
        (
            "head1",
            _to_float_or_none(boltz_results.get("affinity_pred_value1")),
            _to_float_or_none(boltz_results.get("affinity_probability_binary1")),
        ),
        (
            "head2",
            _to_float_or_none(boltz_results.get("affinity_pred_value2")),
            _to_float_or_none(boltz_results.get("affinity_probability_binary2")),
        ),
    ]
    valid_heads = [(name, value, prob) for name, value, prob in heads if value is not None]

    enabled = bool(params.get("affinity_consensus_enabled", False))
    mode = str(params.get("affinity_consensus_mode", "weighted")).strip().lower()
    if mode not in {"weighted", "median", "trimmed_mean"}:
        mode = "weighted"
    weight_floor = max(0.0, float(params.get("affinity_consensus_weight_floor", 0.05)))
    entropy_alpha = min(1.0, max(0.0, float(params.get("affinity_consensus_entropy_alpha", 0.2))))

    if not valid_heads:
        return {
            "selected_affinity_pred_value": 0.0,
            "selected_affinity_probability": 0.0,
            "raw_affinity_pred_value": raw_value,
            "raw_affinity_probability": raw_prob,
            "consensus_enabled": enabled,
            "consensus_mode": mode,
            "consensus_n_heads": 0,
            "consensus_pred_value": None,
            "consensus_probability": None,
            "consensus_std": None,
            "consensus_head_values": {},
            "consensus_head_probabilities": {},
        }

    values = np.array([v for _, v, _ in valid_heads], dtype=float)
    probs = np.array(
        [
            (p if p is not None else 0.5)
            for _, _, p in valid_heads
        ],
        dtype=float,
    )

    if not enabled or len(values) == 1:
        consensus_value = float(values[0] if raw_value is None else raw_value)
        consensus_prob = float(probs[0] if raw_prob is None else raw_prob)
    elif mode == "median":
        consensus_value = float(np.median(values))
        consensus_prob = float(np.median(probs))
    elif mode == "trimmed_mean":
        if len(values) >= 3:
            order = np.argsort(values)
            keep = order[1:-1]
            consensus_value = float(np.mean(values[keep]))
            consensus_prob = float(np.mean(probs[keep]))
        else:
            consensus_value = float(np.mean(values))
            consensus_prob = float(np.mean(probs))
    else:
        # Weighted mode: down-weight uncertain heads where binary probability is near 0.5.
        prob_center_distance = np.abs(probs - 0.5) * 2.0
        entropy = -(
            probs * np.log(np.clip(probs, 1e-8, 1.0))
            + (1.0 - probs) * np.log(np.clip(1.0 - probs, 1e-8, 1.0))
        ) / np.log(2.0)
        weights = np.maximum(weight_floor, prob_center_distance) * (1.0 - entropy_alpha * entropy)
        if float(np.sum(weights)) <= 1e-8:
            weights = np.ones_like(values, dtype=float)
        consensus_value = float(np.average(values, weights=weights))
        consensus_prob = float(np.average(probs, weights=weights))

    consensus_std = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    selected_value = consensus_value if (enabled and consensus_value is not None) else raw_value
    selected_prob = consensus_prob if (enabled and consensus_prob is not None) else raw_prob
    if selected_value is None:
        selected_value = float(values[0])
    if selected_prob is None:
        selected_prob = float(probs[0])

    return {
        "selected_affinity_pred_value": float(selected_value),
        "selected_affinity_probability": float(selected_prob),
        "raw_affinity_pred_value": raw_value,
        "raw_affinity_probability": raw_prob,
        "consensus_enabled": enabled,
        "consensus_mode": mode,
        "consensus_n_heads": int(len(values)),
        "consensus_pred_value": float(consensus_value) if consensus_value is not None else None,
        "consensus_probability": float(consensus_prob) if consensus_prob is not None else None,
        "consensus_std": float(consensus_std) if consensus_std is not None else None,
        "consensus_head_values": {name: value for name, value, _ in valid_heads},
        "consensus_head_probabilities": {
            name: prob for name, _, prob in valid_heads if prob is not None
        },
    }


def _create_result_entry(
    protein_name: str,
    protein_sequence: str,
    drug_name: str,
    smiles: str,
    workspace_name: str,
    design_name: str,
    boltz_results: Dict[str, Any],
    params: Dict[str, Any],
) -> Dict[str, Any]:
    affinity_meta = _compute_affinity_consensus(boltz_results=boltz_results, params=params)
    affinity_pred_value = affinity_meta["selected_affinity_pred_value"]
    affinity_probability_binary = affinity_meta["selected_affinity_probability"]
    confidence = boltz_results.get("confidence_score", 0.85)
    ptm = boltz_results.get("ptm", 0.8)
    iptm = boltz_results.get("iptm", 0.7)
    avg_plddt = boltz_results.get("complex_plddt", 0.85) * 100

    ic50 = 10 ** (affinity_pred_value)
    pic50 = -np.log10(ic50 * 1e-6) if ic50 > 0 else None

    return {
        "protein_name": protein_name,
        "drug_name": drug_name,
        "protein_sequence": protein_sequence,
        "smiles": smiles,
        "ic50_um": ic50,
        "pic50": pic50,
        "affinity_probability": affinity_probability_binary,
        "affinity_pred_value_raw": affinity_meta["raw_affinity_pred_value"],
        "affinity_probability_raw": affinity_meta["raw_affinity_probability"],
        "affinity_pred_value_consensus": affinity_meta["consensus_pred_value"],
        "affinity_probability_consensus": affinity_meta["consensus_probability"],
        "affinity_consensus_std": affinity_meta["consensus_std"],
        "affinity_consensus_n_heads": affinity_meta["consensus_n_heads"],
        "affinity_consensus_mode": affinity_meta["consensus_mode"],
        "affinity_consensus_enabled": affinity_meta["consensus_enabled"],
        "affinity_pred_value_multisampling_aggregate": boltz_results.get("affinity_multisampling_aggregate"),
        "affinity_probability_multisampling_aggregate": boltz_results.get("affinity_multisampling_probability_aggregate"),
        "affinity_pred_value1": boltz_results.get("affinity_pred_value1"),
        "affinity_probability1": boltz_results.get("affinity_probability_binary1"),
        "affinity_pred_value2": boltz_results.get("affinity_pred_value2"),
        "affinity_probability2": boltz_results.get("affinity_probability_binary2"),
        "affinity_multisampling_applied": boltz_results.get("affinity_multisampling_applied"),
        "affinity_multisampling_mode": boltz_results.get("affinity_multisampling_mode"),
        "affinity_multisampling_n": boltz_results.get("affinity_multisampling_n"),
        "affinity_multisampling_std": boltz_results.get("affinity_multisampling_std"),
        "affinity_multisampling_min": boltz_results.get("affinity_multisampling_min"),
        "affinity_multisampling_max": boltz_results.get("affinity_multisampling_max"),
        "affinity_multisampling_setting_values": boltz_results.get("affinity_multisampling_setting_values"),
        "affinity_multisampling_setting_probabilities": boltz_results.get("affinity_multisampling_setting_probabilities"),
        "affinity_multisampling_refinement_steps": boltz_results.get("affinity_multisampling_refinement_steps"),
        "affinity_multisampling_diffusion_samples": boltz_results.get("affinity_multisampling_diffusion_samples"),
        "affinity_multisampling_profiles": boltz_results.get("affinity_multisampling_setting_profiles"),
        "affinity_multisampling_setting_profiles": boltz_results.get("affinity_multisampling_setting_profiles"),
        "affinity_multisampling_recommended_setting": boltz_results.get("affinity_multisampling_recommended_setting"),
        "affinity_multisampling_weighting_mode": boltz_results.get("affinity_multisampling_weighting_mode"),
        "affinity_multisampling_setting_uncertainties": boltz_results.get("affinity_multisampling_setting_uncertainties"),
        "affinity_multisampling_setting_weights": boltz_results.get("affinity_multisampling_setting_weights"),
        "affinity_patch_mode": boltz_results.get("affinity_patch_mode"),
        "affinity_patch_applied": boltz_results.get("affinity_patch_applied"),
        "affinity_patch_topk_applied": boltz_results.get("affinity_patch_topk_applied"),
        "affinity_patch_topk_k": boltz_results.get("affinity_patch_topk_k"),
        "affinity_patch_best_idx_raw": boltz_results.get("affinity_patch_best_idx_raw"),
        "affinity_patch_topk_indices": boltz_results.get("affinity_patch_topk_indices"),
        "affinity_patch_topk_weights": boltz_results.get("affinity_patch_topk_weights"),
        "affinity_patch_selector_scores": boltz_results.get("affinity_patch_selector_scores"),
        "affinity_patch_selector_local_scores": boltz_results.get("affinity_patch_selector_local_scores"),
        "affinity_patch_selector_local_plddt": boltz_results.get("affinity_patch_selector_local_plddt"),
        "affinity_patch_selector_local_ipde": boltz_results.get("affinity_patch_selector_local_ipde"),
        "affinity_patch_selector_mutation_positions": boltz_results.get("affinity_patch_selector_mutation_positions"),
        "confidence": confidence,
        "ptm": ptm,
        "iptm": iptm,
        "avg_plddt": avg_plddt,
        "status": "Success",
        "workspace": workspace_name,
        "design": design_name,
        "cofactor_info": params.get("cofactor_info"),
        "boltz2_parameters": {
            "use_gpu": params.get("use_gpu", True),
            "accelerator": params.get("accelerator", "gpu"),
            "devices": params.get("devices", 1),
            "cuda_visible_devices": params.get("cuda_visible_devices"),
            "preprocessing_threads": params.get("preprocessing_threads", 1),
            "recycling_steps": params.get("recycling_steps", 3),
            "sampling_steps": params.get("sampling_steps", 200),
            "diffusion_samples": params.get("diffusion_samples", 1),
            "max_parallel_samples": params.get("max_parallel_samples", 5),
            "step_scale": params.get("step_scale", 1.638),
            "affinity_mw_correction": params.get("affinity_mw_correction", False),
            "affinity_consensus_enabled": params.get("affinity_consensus_enabled", False),
            "affinity_consensus_mode": params.get("affinity_consensus_mode", "weighted"),
            "affinity_consensus_weight_floor": params.get("affinity_consensus_weight_floor", 0.05),
            "affinity_consensus_entropy_alpha": params.get("affinity_consensus_entropy_alpha", 0.2),
            "affinity_multisampling_enabled": params.get("affinity_multisampling_enabled", False),
            "affinity_multisampling_profiles": params.get("affinity_multisampling_profiles"),
            "affinity_multisampling_settings": params.get("affinity_multisampling_settings"),
            "affinity_multisampling_refinement_steps": params.get("affinity_multisampling_refinement_steps"),
            "affinity_multisampling_aggregate_mode": params.get("affinity_multisampling_aggregate_mode", "median"),
            "affinity_multisampling_apply_aggregate": params.get("affinity_multisampling_apply_aggregate", True),
            "affinity_multisampling_early_stop_enabled": params.get("affinity_multisampling_early_stop_enabled", True),
            "affinity_multisampling_early_stop_min_points": params.get("affinity_multisampling_early_stop_min_points", 2),
            "affinity_multisampling_early_stop_delta": params.get("affinity_multisampling_early_stop_delta", 0.02),
            "affinity_multisampling_early_stop_std": params.get("affinity_multisampling_early_stop_std", 0.04),
            "affinity_multisampling_early_stop_patience": params.get("affinity_multisampling_early_stop_patience", 1),
            "affinity_multisampling_robust_outlier_filter": params.get("affinity_multisampling_robust_outlier_filter", True),
            "affinity_multisampling_robust_outlier_zmax": params.get("affinity_multisampling_robust_outlier_zmax", 3.5),
            "affinity_multisampling_bootstrap_samples": params.get("affinity_multisampling_bootstrap_samples", 300),
            "external_boltz_patch_enabled": params.get("external_boltz_patch_enabled", False),
            "external_boltz_patch_mode": params.get("external_boltz_patch_mode", "mutation_aware_v2"),
            "external_boltz_patch_weight_floor": params.get("external_boltz_patch_weight_floor", 0.05),
            "external_boltz_patch_entropy_alpha": params.get("external_boltz_patch_entropy_alpha", 0.20),
            "external_boltz_patch_uncertainty_penalty": params.get("external_boltz_patch_uncertainty_penalty", 0.15),
            "external_boltz_patch_min_confidence": params.get("external_boltz_patch_min_confidence", 0.35),
            "use_potentials": params.get("use_potentials", False),
            "method": params.get("method"),
            "max_msa_seqs": params.get("max_msa_seqs", 8192),
            "sampling_steps_affinity": params.get("sampling_steps_affinity", 200),
            "diffusion_samples_affinity": params.get("diffusion_samples_affinity", 5),
            "subsample_msa": params.get("subsample_msa", False),
            "num_subsampled_msa": params.get("num_subsampled_msa", 1024),
            "enable_retries": params.get("enable_retries", True),
            "max_retry_attempts": params.get("max_retry_attempts", 2),
            "retry_delay_base": params.get("retry_delay_base", 5),
            "template_cif_path": params.get("template_cif_path"),
            "mutation_steering_config": params.get("mutation_steering_config"),
            "use_cached_msa": params.get("use_cached_msa", False),
            "enable_msa_cache": params.get("enable_msa_cache", True),
            "override": params.get("override", False),
            "prediction_timeout_seconds": params.get("prediction_timeout_seconds", 300),
        },
    }


def _persist_result_entry(
    project_name: str,
    result_entry: Dict[str, Any],
    computation_time: Optional[float],
    binding_pocket_constraints: Optional[Dict[str, Any]] = None,
    template_cif_path: Optional[str] = None,
    boltz_command: Optional[str] = None,
    log_warning: bool = True,
) -> None:
    if not project_name or not result_entry:
        return
    boltz_commands = [boltz_command] if boltz_command else None
    try:
        save_screening_results(
            [result_entry],
            project_name,
            RESULTS_DIR,
            computation_time=computation_time,
            template_cif_path=template_cif_path,
            binding_pocket_constraints=binding_pocket_constraints,
            boltz_commands=boltz_commands,
        )
    except Exception as exc:
        message = f"Unable to persist screening result: {str(exc)[:150]}"
        if log_warning:
            st.warning(f":material/warning: {message}")
        else:
            print(f"[WARN] {message}")


def execute_screening_job(job: ScreeningJob, worker_id: int = 0) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    params = job.parameters or {}
    ligand_smiles = "" if job.structure_only else job.smiles
    ligand_display_name = "" if job.structure_only else job.drug_name
    start_time = time.time()
    enable_msa_cache = bool(params.get("enable_msa_cache", True))
    wt_sequence_for_msa_cache = params.get("wt_sequence_for_msa_cache")
    msa_path, use_cached_msa = resolve_msa_cache_for_prediction(
        protein_sequence=job.protein_sequence,
        wt_sequence=wt_sequence_for_msa_cache,
        enable_msa_cache=enable_msa_cache,
    )
    params["use_cached_msa"] = use_cached_msa
    if msa_path:
        params["msa_path"] = msa_path

    queue_gpu_devices = params.get("queue_gpu_devices")
    effective_cuda_visible_devices = params.get("cuda_visible_devices")
    if (
        params.get("accelerator", "gpu") == "gpu"
        and isinstance(queue_gpu_devices, list)
        and len(queue_gpu_devices) > 0
    ):
        selected_idx = int(worker_id) % len(queue_gpu_devices)
        effective_cuda_visible_devices = str(queue_gpu_devices[selected_idx])
        params["cuda_visible_devices"] = effective_cuda_visible_devices
        params["devices"] = 1

    boltz_results, success, error_message = run_boltz_with_retry(
        workspace_name=job.workspace_name,
        design_name=job.design_name,
        protein_sequence=job.protein_sequence,
        ligand_smiles=ligand_smiles,
        project_name=job.project_name,
        protein_display_name=job.protein_name,
        ligand_display_name=ligand_display_name,
        use_gpu=params.get("use_gpu", True),
        binding_pocket_constraints=params.get("binding_pocket_constraints"),
        override=params.get("override", False),
        recycling_steps=params.get("recycling_steps", 3),
        sampling_steps=params.get("sampling_steps", 200),
        diffusion_samples=params.get("diffusion_samples", 1),
        max_parallel_samples=params.get("max_parallel_samples", 5),
        step_scale=params.get("step_scale", 1.638),
        affinity_mw_correction=params.get("affinity_mw_correction", False),
        external_boltz_patch_enabled=params.get("external_boltz_patch_enabled", False),
        external_boltz_patch_mode=params.get("external_boltz_patch_mode", "mutation_aware_v2"),
        external_boltz_patch_weight_floor=params.get("external_boltz_patch_weight_floor", 0.05),
        external_boltz_patch_entropy_alpha=params.get("external_boltz_patch_entropy_alpha", 0.20),
        external_boltz_patch_uncertainty_penalty=params.get("external_boltz_patch_uncertainty_penalty", 0.15),
        external_boltz_patch_min_confidence=params.get("external_boltz_patch_min_confidence", 0.35),
        affinity_multisampling_enabled=params.get("affinity_multisampling_enabled", False),
        affinity_multisampling_profiles=params.get("affinity_multisampling_profiles"),
        affinity_multisampling_settings=params.get("affinity_multisampling_settings"),
        affinity_multisampling_refinement_steps=params.get("affinity_multisampling_refinement_steps"),
        affinity_multisampling_aggregate_mode=params.get("affinity_multisampling_aggregate_mode", "median"),
        affinity_multisampling_apply_aggregate=params.get("affinity_multisampling_apply_aggregate", True),
        affinity_multisampling_early_stop_enabled=params.get("affinity_multisampling_early_stop_enabled", True),
        affinity_multisampling_early_stop_min_points=params.get("affinity_multisampling_early_stop_min_points", 2),
        affinity_multisampling_early_stop_delta=params.get("affinity_multisampling_early_stop_delta", 0.02),
        affinity_multisampling_early_stop_std=params.get("affinity_multisampling_early_stop_std", 0.04),
        affinity_multisampling_early_stop_patience=params.get("affinity_multisampling_early_stop_patience", 1),
        affinity_multisampling_robust_outlier_filter=params.get("affinity_multisampling_robust_outlier_filter", True),
        affinity_multisampling_robust_outlier_zmax=params.get("affinity_multisampling_robust_outlier_zmax", 3.5),
        affinity_multisampling_bootstrap_samples=params.get("affinity_multisampling_bootstrap_samples", 300),
        confidence_target=params.get("confidence_target", "balanced"),
        max_msa_seqs=params.get("max_msa_seqs", 8192),
        sampling_steps_affinity=params.get("sampling_steps_affinity", 200),
        diffusion_samples_affinity=params.get("diffusion_samples_affinity", 5),
        cofactor_info=params.get("cofactor_info"),
        enable_retries=params.get("enable_retries", True),
        max_retry_attempts=params.get("max_retry_attempts", 2),
        retry_delay_base=params.get("retry_delay_base", 5),
        subsample_msa=params.get("subsample_msa", False),
        num_subsampled_msa=params.get("num_subsampled_msa", 1024),
        template_cif_path=params.get("template_cif_path"),
        structure_only=job.structure_only,
        ptm_modifications=params.get("ptm_modifications"),
        prediction_timeout_seconds=params.get("prediction_timeout_seconds", 300),
        accelerator=params.get("accelerator", "gpu"),
        devices=params.get("devices", 1),
        cuda_visible_devices=effective_cuda_visible_devices,
        preprocessing_threads=params.get("preprocessing_threads", 1),
        use_potentials=params.get("use_potentials", False),
        method=params.get("method"),
        emit_streamlit_feedback=False,
        msa_path=msa_path,
        use_cached_msa=use_cached_msa,
        enable_msa_cache=enable_msa_cache,
    )
    if not success or boltz_results is None:
        raise RuntimeError(error_message or "Boltz prediction failed")

    result_entry = _create_result_entry(
        protein_name=job.protein_name,
        protein_sequence=job.protein_sequence,
        drug_name=job.drug_name if not job.structure_only else "",
        smiles=job.smiles if not job.structure_only else "",
        workspace_name=job.workspace_name,
        design_name=job.design_name,
        boltz_results=boltz_results,
        params=params,
    )

    computation_time = time.time() - start_time
    timestamp = datetime.now().isoformat()
    metadata = {
        "computation_time_seconds": computation_time,
        "timestamp": timestamp,
        "worker_id": worker_id,
        "worker_cuda_visible_devices": effective_cuda_visible_devices,
        "binding_pocket_constraints": params.get("binding_pocket_constraints"),
        "template_cif_path": params.get("template_cif_path"),
        "boltz_command": _build_boltz_command(yaml_filename, params),
    }

    _persist_result_entry(
        job.project_name,
        result_entry,
        computation_time,
        binding_pocket_constraints=params.get("binding_pocket_constraints"),
        template_cif_path=params.get("template_cif_path"),
        boltz_command=metadata["boltz_command"],
        log_warning=False,
    )
    metadata["persisted"] = True
    return result_entry, metadata


def ensure_job_manager_executor(worker_count: int = 1) -> None:
    if not USE_SCREENING_JOB_QUEUE:
        return
    manager = get_job_manager()
    if manager is None:
        return
    session_state = getattr(st, "session_state", None)
    already_registered = bool(session_state is not None and session_state.get(JOB_MANAGER_EXECUTOR_STATE_KEY))
    if not already_registered:
        manager.register_executor(execute_screening_job)
        if session_state is not None:
            session_state[JOB_MANAGER_EXECUTOR_STATE_KEY] = True
    manager.set_worker_count(max(1, int(worker_count)))


def prepare_screening_jobs(
    protein_sequences: List[Tuple[str, str]],
    drug_smiles: List[Tuple[str, str]],
    project_name: str,
    structure_only: bool,
    use_existing_results: bool,
    protein_drug_filter: Optional[Dict[str, Any]],
    shared_params: Dict[str, Any],
    manager: Optional[ScreeningJobManager] = None,
) -> Tuple[List[ScreeningJob], List[Dict[str, Any]], Dict[str, Any]]:
    if manager is None and USE_SCREENING_JOB_QUEUE:
        manager = get_job_manager()
    jobs: List[ScreeningJob] = []
    cached_results: List[Dict[str, Any]] = []
    summary = {
        "existing_results": 0,
        "new_jobs": 0,
        "skipped": 0,
        "duplicate_jobs": 0,
        "warnings": [],
    }
    seen_signatures: Set[str] = set()
    wt_sequence_for_msa_cache = next(
        (seq for name, seq in protein_sequences if str(name).strip().upper() == "WT"),
        None
    )

    ordered_protein_sequences = _order_protein_sequences_for_screening(protein_sequences)

    for protein_name, protein_seq in ordered_protein_sequences:
        if structure_only:
            if not should_evaluate_protein_drug_pair(protein_name, None, protein_drug_filter):
                summary["skipped"] += 1
                continue
            combos = [("", "")]
        else:
            combos = drug_smiles

        for drug_name, drug_smiles_str in combos:
            if not structure_only and not should_evaluate_protein_drug_pair(protein_name, drug_name, protein_drug_filter):
                summary["skipped"] += 1
                continue

            params = copy.deepcopy(shared_params)
            params["cofactor_info"] = shared_params.get("cofactor_info")
            params["mutation_steering_config"] = shared_params.get("mutation_steering_config")
            params["binding_pocket_constraints"] = _resolve_mutation_steering_constraints(
                protein_name=protein_name,
                protein_sequence=protein_seq,
                base_constraints=shared_params.get("binding_pocket_constraints"),
                structure_only=structure_only,
                mutation_steering_config=params.get("mutation_steering_config"),
            )
            params["ptm_modifications"] = shared_params.get("ptm_modifications")
            params["queue_gpu_devices"] = shared_params.get("queue_gpu_devices")
            params["enable_msa_cache"] = bool(shared_params.get("enable_msa_cache", True))
            params["structure_only"] = structure_only
            if wt_sequence_for_msa_cache:
                params["wt_sequence_for_msa_cache"] = wt_sequence_for_msa_cache

            workspace_name, design_name, stable_job_id = _derive_stable_job_identifiers(
                project_name=project_name,
                protein_name=protein_name,
                protein_sequence=protein_seq,
                drug_name="" if structure_only else drug_name,
                smiles="" if structure_only else drug_smiles_str,
                structure_only=structure_only,
                params=params,
            )

            if use_existing_results:
                existing_yaml_name = _find_existing_results_by_yaml_name(
                    project_name=project_name,
                    yaml_name=f"{workspace_name}_{design_name}",
                    structure_only=structure_only,
                )
                if not existing_yaml_name:
                    existing_yaml_name = find_existing_screening_results(protein_name, drug_name, project_name)
                if existing_yaml_name:
                    project_dir = os.path.join(RESULTS_DIR, project_name)
                    yaml_filepath = os.path.join(project_dir, f"{existing_yaml_name}.yaml")
                    is_valid, validation_error = validate_boltz_results(yaml_filepath, structure_only=structure_only)
                    if is_valid:
                        boltz_results = get_screening_existing_boltz_results(
                            existing_yaml_name,
                            project_name,
                            binding_pocket_constraints=params.get("binding_pocket_constraints"),
                            cofactor_info=params.get("cofactor_info"),
                            structure_only=structure_only,
                        )
                        if boltz_results:
                            workspace_name, design_name = existing_yaml_name, existing_yaml_name
                            result_entry = _create_result_entry(
                                protein_name=protein_name,
                                protein_sequence=protein_seq,
                                drug_name=drug_name if not structure_only else "",
                                smiles=drug_smiles_str if not structure_only else "",
                                workspace_name=workspace_name,
                                design_name=design_name,
                                boltz_results=boltz_results,
                                params=params,
                            )
                            cached_results.append(result_entry)
                            summary["existing_results"] += 1
                            continue
                        else:
                            summary["warnings"].append(
                                f"Could not parse existing results for {protein_name}{(' + ' + drug_name) if drug_name else ''}"
                            )
                    else:
                        summary["warnings"].append(
                            f"Existing results invalid for {protein_name}{(' + ' + drug_name) if drug_name else ''}: {validation_error}"
                        )

            job_params = copy.deepcopy(params)
            max_attempts = int(job_params.get("max_retry_attempts", 0)) + 1 if job_params.get("enable_retries", True) else 1

            normalized_drug_name = drug_name if not structure_only else ""
            normalized_smiles = drug_smiles_str if not structure_only else ""
            signature = ScreeningJob.build_signature(
                project_name=project_name,
                protein_name=protein_name,
                protein_sequence=protein_seq,
                drug_name=normalized_drug_name,
                smiles=normalized_smiles,
                structure_only=structure_only,
                parameters=job_params,
            )
            ligand_display = normalized_drug_name or ("No ligand" if structure_only else "")
            combo_label = f"{protein_name}{' + ' + ligand_display if ligand_display else ''}"

            if signature in seen_signatures:
                summary["duplicate_jobs"] += 1
                summary["warnings"].append(f"Skipped duplicate configuration requested for {combo_label}.")
                continue

            existing_job = manager.has_job_with_signature(signature) if manager else None
            if existing_job:
                summary["duplicate_jobs"] += 1
                summary["warnings"].append(
                    f"Skipped {combo_label} because an equivalent job already exists (status: {existing_job.status})."
                )
                continue

            seen_signatures.add(signature)
            job_id = stable_job_id
            job = ScreeningJob(
                job_id=job_id,
                project_name=project_name,
                protein_name=protein_name,
                protein_sequence=protein_seq,
                drug_name=normalized_drug_name,
                smiles=normalized_smiles,
                structure_only=structure_only,
                parameters=job_params,
                workspace_name=workspace_name,
                design_name=design_name,
                max_attempts=max_attempts,
            )
            job.set_cached_signature(signature)
            jobs.append(job)
            summary["new_jobs"] += 1

    return jobs, cached_results, summary


def synchronize_job_results(project_name: str) -> None:
    if not USE_SCREENING_JOB_QUEUE or not project_name:
        return
    manager = get_job_manager()
    if manager is None:
        return
    ready = manager.get_uncommitted_results(project_name)
    if not ready:
        return

    new_results: List[Dict[str, Any]] = []
    job_ids: List[str] = []
    computation_times: List[float] = []

    for job, result, metadata in ready:
        new_results.append(result)
        job_ids.append(job.job_id)
        duration = metadata.get("computation_time_seconds") if metadata else None
        if isinstance(duration, (int, float)):
            computation_times.append(duration)
        if not metadata or not metadata.get("persisted"):
            _persist_result_entry(
                project_name,
                result,
                duration,
                binding_pocket_constraints=metadata.get("binding_pocket_constraints") if metadata else None,
                template_cif_path=metadata.get("template_cif_path") if metadata else None,
                boltz_command=metadata.get("boltz_command") if metadata else None,
            )

    if new_results:
        current_results = st.session_state.get("screening_results", [])
        combined_results = deduplicate_results(current_results + new_results)
        st.session_state.screening_results = combined_results
        if computation_times:
            st.session_state.last_computation_time = sum(computation_times) / len(computation_times)

    manager.mark_results_committed(job_ids)


def render_job_queue_status(project_name: str) -> Optional[Dict[str, Any]]:
    if not USE_SCREENING_JOB_QUEUE or not project_name:
        return None
    manager = get_job_manager()
    if manager is None:
        try:
            session_state = st.session_state
        except Exception:
            session_state = None
        if session_state is not None:
            init_error = session_state.get(JOB_MANAGER_INIT_ERROR_KEY)
        else:
            init_error = None
        if init_error:
            st.error(f"Job queue unavailable: {init_error}")
        return None
    summary = manager.get_project_summary(project_name)
    total = summary.get("total", 0)
    if total == 0:
        return summary
    completed = summary.get("success", 0) + summary.get("failed", 0) + summary.get("cancelled", 0)
    progress = completed / total if total else 0

    completed_display = completed
    status_text = (
        f"{summary.get('running', 0)} running • "
        f"{summary.get('pending', 0)} pending • "
        f"{completed_display} complete"
    )
    if summary.get("failed", 0):
        status_text += f" ({summary.get('failed', 0)} failed)"
    st.progress(progress, text=status_text)

    elapsed_display = _format_duration(summary.get("elapsed_seconds"))
    eta_display = _format_duration(summary.get("eta_seconds"))
    timing_messages = []
    if elapsed_display:
        timing_messages.append(f"Elapsed {elapsed_display}")
    if eta_display:
        timing_messages.append(f"ETA {eta_display}")
    if timing_messages:
        st.caption(" • ".join(timing_messages))

    active_jobs = summary.get("active_jobs") or []
    if active_jobs:
        labels: List[str] = []
        for job in active_jobs[:3]:
            ligand_display = job.drug_name or ("No ligand" if job.structure_only else "")
            label = f"{job.protein_name}{' + ' + ligand_display if ligand_display else ''}"
            labels.append(label)
        if len(active_jobs) > 3:
            labels.append(f"... and {len(active_jobs) - 3} more")
        st.info(f"Processing: {' | '.join(labels)}")

    project_jobs = manager.get_project_jobs(project_name)
    if project_jobs:
        job_rows: List[Dict[str, Any]] = []
        for job in project_jobs:
            ligand_display = job.drug_name or ("No ligand" if job.structure_only else "")
            started_at = datetime.fromtimestamp(job.started_at).strftime("%Y-%m-%d %H:%M:%S") if job.started_at else ""
            completed_at = datetime.fromtimestamp(job.completed_at).strftime("%Y-%m-%d %H:%M:%S") if job.completed_at else ""
            created_at = datetime.fromtimestamp(job.created_at).strftime("%Y-%m-%d %H:%M:%S") if job.created_at else ""
            duration = ""
            if job.started_at and job.completed_at and job.completed_at >= job.started_at:
                duration_seconds = job.completed_at - job.started_at
                duration = _format_duration(duration_seconds)
            job_rows.append(
                {
                    "Protein": job.protein_name,
                    "Ligand": ligand_display,
                    "Status": job.status,
                    "Retries": job.retries,
                    "Created": created_at,
                    "Started": started_at,
                    "Completed": completed_at,
                    "Elapsed": duration,
                    "Error": job.error or "",
                }
            )
        with st.expander("Job Queue Details", expanded=False, icon=":material/list:"):
            st.dataframe(
                pd.DataFrame(job_rows),
                use_container_width=True,
                hide_index=True,
            )
    return summary


def maybe_schedule_queue_autorefresh(
    summary: Optional[Dict[str, Any]],
    interval_seconds: float = QUEUE_STATUS_REFRESH_INTERVAL_SECONDS,
) -> None:
    if not summary:
        try:
            state = st.session_state
        except Exception:
            state = None
        if state is not None and QUEUE_NEXT_REFRESH_STATE_KEY in state:
            del state[QUEUE_NEXT_REFRESH_STATE_KEY]
        return
    if interval_seconds is None or interval_seconds <= 0:
        return
    has_active_jobs = (summary.get("running", 0) or 0) > 0 or (summary.get("pending", 0) or 0) > 0
    try:
        state = st.session_state
    except Exception:
        state = None
    if not has_active_jobs:
        if state is not None and QUEUE_NEXT_REFRESH_STATE_KEY in state:
            del state[QUEUE_NEXT_REFRESH_STATE_KEY]
        return
    if state is None:
        return
    now = time.time()
    next_refresh_at = state.get(QUEUE_NEXT_REFRESH_STATE_KEY, 0.0)
    if now < next_refresh_at:
        return
    state[QUEUE_NEXT_REFRESH_STATE_KEY] = now + interval_seconds
    time.sleep(interval_seconds)
    try:
        st.rerun()
    except Exception:
        pass


def check_screening_existing_boltz_results(workspace_name, design_name, project_name):
    """Check if existing Boltz results are available for the given workspace and design name in project-specific directory."""
    try:
        # Create the expected YAML filename
        yaml_filename = f"{workspace_name}_{design_name}.yaml"
        yaml_name = os.path.splitext(yaml_filename)[0]
        
        # Check if project directory exists
        project_dir = os.path.join("boltzomics_screening_results", project_name)
        if not os.path.exists(project_dir):
            return False
        
        # Check if the specific results folder exists
        boltz_results_dir = os.path.join(project_dir, f"boltz_results_{yaml_name}")
        if not os.path.exists(boltz_results_dir):
            return False
        
        # Check if both affinity and confidence files exist
        affinity_file = os.path.join(boltz_results_dir, "predictions", yaml_name, f"affinity_{yaml_name}.json")
        confidence_file = os.path.join(boltz_results_dir, "predictions", yaml_name, f"confidence_{yaml_name}_model_0.json")
        
        # Check if both files exist and are not empty
        if not os.path.exists(affinity_file) or not os.path.exists(confidence_file):
            return False
        
        # Check if files have content (not empty)
        if os.path.getsize(affinity_file) == 0 or os.path.getsize(confidence_file) == 0:
            return False
        
        return True
    except Exception:
        return False

def find_existing_screening_results_manual(protein_name, drug_name, project_name):
    """
    Manual method to find existing screening results by comprehensive directory and file checking.
    This method provides more reliable detection than boltz2-based checking.
    """
    try:
        project_dir = os.path.join("boltzomics_screening_results", project_name)
        if not os.path.exists(project_dir):
            return None
        
        # Create the design name pattern
        design_name = f"{protein_name}_{drug_name}".replace(' ', '_').replace(':', '_')
        
        # First check for exact YAML files matching the pattern
        for item in os.listdir(project_dir):
            if item.endswith('.yaml') and design_name in item:
                yaml_name = item[:-5]  # Remove .yaml extension
                yaml_path = os.path.join(project_dir, item)
                
                # Verify YAML file is valid and not empty
                if os.path.exists(yaml_path) and os.path.getsize(yaml_path) > 0:
                    # Check for corresponding boltz_results directory
                    results_dir = os.path.join(project_dir, f"boltz_results_{yaml_name}")
                    if os.path.exists(results_dir):
                        # Verify required prediction files exist
                        predictions_dir = os.path.join(results_dir, "predictions", yaml_name)
                        if os.path.exists(predictions_dir):
                            affinity_file = os.path.join(predictions_dir, f"affinity_{yaml_name}.json")
                            confidence_file = os.path.join(predictions_dir, f"confidence_{yaml_name}_model_0.json")
                            
                            # Check if at least one of the required files exists and is not empty
                            affinity_exists = os.path.exists(affinity_file) and os.path.getsize(affinity_file) > 0
                            confidence_exists = os.path.exists(confidence_file) and os.path.getsize(confidence_file) > 0
                            
                            if affinity_exists or confidence_exists:
                                return yaml_name
        
        # If no exact match found, check for partial matches in directory names
        # This handles cases where directory names might have slight variations
        for item in os.listdir(project_dir):
            if item.startswith("boltz_results_screening_"):
                # Extract components and check if they match
                dir_parts = item.replace("boltz_results_screening_", "").split("_")
                if len(dir_parts) >= 3:  # timestamp_protein_drug format expected
                    # Reconstruct protein and drug names from directory
                    potential_protein = "_".join(dir_parts[1:-1]) if len(dir_parts) > 3 else dir_parts[1]
                    potential_drug = dir_parts[-1]
                    
                    # Check if protein and drug names match (case-insensitive, flexible matching)
                    protein_clean = protein_name.replace(' ', '_').replace(':', '_').lower()
                    drug_clean = drug_name.replace(' ', '_').replace(':', '_').lower()
                    
                    if (potential_protein.lower() == protein_clean and 
                        potential_drug.lower() == drug_clean):
                        
                        yaml_name = item.replace("boltz_results_", "")
                        yaml_path = os.path.join(project_dir, f"{yaml_name}.yaml")
                        
                        # Verify YAML file exists
                        if os.path.exists(yaml_path) and os.path.getsize(yaml_path) > 0:
                            # Check for required prediction files
                            predictions_dir = os.path.join(project_dir, item, "predictions", yaml_name)
                            if os.path.exists(predictions_dir):
                                affinity_file = os.path.join(predictions_dir, f"affinity_{yaml_name}.json")
                                confidence_file = os.path.join(predictions_dir, f"confidence_{yaml_name}_model_0.json")
                                
                                affinity_exists = os.path.exists(affinity_file) and os.path.getsize(affinity_file) > 0
                                confidence_exists = os.path.exists(confidence_file) and os.path.getsize(confidence_file) > 0
                                
                                if affinity_exists or confidence_exists:
                                    return yaml_name
        
        return None
    except Exception as e:
        # Log the exception for debugging but don't crash
        return None

def find_existing_screening_results(protein_name, drug_name, project_name):
    """
    Find existing screening results by trying both manual checking and boltz2-based checking.
    Manual checking is performed first as it's more reliable.
    """
    try:
        # First try manual checking method
        manual_result = find_existing_screening_results_manual(protein_name, drug_name, project_name)
        if manual_result:
            return manual_result
        
        # Fallback to original boltz2-based method for backward compatibility
        project_dir = os.path.join("boltzomics_screening_results", project_name)
        if not os.path.exists(project_dir):
            return None
        
        # Create the design name pattern
        design_name = f"{protein_name}_{drug_name}".replace(' ', '_').replace(':', '_')
        
        # Search for existing boltz_results directories that match the design name
        for item in os.listdir(project_dir):
            if item.startswith("boltz_results_screening_") and design_name in item:
                # Extract the full yaml name from the directory
                yaml_name = item.replace("boltz_results_", "")
                
                # Check if the required files exist
                affinity_file = os.path.join(project_dir, item, "predictions", yaml_name, f"affinity_{yaml_name}.json")
                confidence_file = os.path.join(project_dir, item, "predictions", yaml_name, f"confidence_{yaml_name}_model_0.json")
                
                if (os.path.exists(affinity_file) and os.path.exists(confidence_file) and 
                    os.path.getsize(affinity_file) > 0 and os.path.getsize(confidence_file) > 0):
                    return yaml_name
        
        return None
    except Exception:
        return None

def get_screening_existing_boltz_results(yaml_name, project_name, binding_pocket_constraints=None, cofactor_info=None, structure_only=False):
    """Get existing Boltz results if available in project-specific directory."""
    try:
        # Create the YAML filepath using the found yaml name
        project_dir = os.path.join("boltzomics_screening_results", project_name)
        yaml_filepath = os.path.join(project_dir, f"{yaml_name}.yaml")
        
        # Check if the YAML file exists
        if not os.path.exists(yaml_filepath):
            return None
        
        # Parse results using the existing YAML file
        return utils.parse_boltz_results(yaml_filepath, structure_only=structure_only)
    except Exception as e:
        st.warning(f"Error getting existing Boltz results: {e}")
        return None

def run_screening_prediction(
    protein_sequences: List[Tuple[str, str]],
    drug_smiles: List[Tuple[str, str]],
    project_name: str,
    use_gpu: bool = True,
    use_existing_results: bool = True,
    recycling_steps: int = 4,
    sampling_steps: int = 300,
    diffusion_samples: int = 1,
    max_parallel_samples: int = 5,
    step_scale: float = 1.638,
    affinity_mw_correction: bool = False,
    affinity_consensus_enabled: bool = False,
    affinity_consensus_mode: str = "weighted",
    affinity_consensus_weight_floor: float = 0.05,
    affinity_consensus_entropy_alpha: float = 0.2,
    external_boltz_patch_enabled: bool = False,
    external_boltz_patch_mode: str = "mutation_aware_v2",
    external_boltz_patch_weight_floor: float = 0.05,
    external_boltz_patch_entropy_alpha: float = 0.20,
    external_boltz_patch_uncertainty_penalty: float = 0.15,
    external_boltz_patch_min_confidence: float = 0.35,
    affinity_multisampling_enabled: bool = False,
    affinity_multisampling_profiles: Optional[List[Dict[str, int]]] = None,
    affinity_multisampling_settings: Optional[List[str]] = None,
    affinity_multisampling_refinement_steps: Optional[List[int]] = None,
    affinity_multisampling_aggregate_mode: str = "median",
    affinity_multisampling_apply_aggregate: bool = True,
    affinity_multisampling_early_stop_enabled: bool = True,
    affinity_multisampling_early_stop_min_points: int = 2,
    affinity_multisampling_early_stop_delta: float = 0.02,
    affinity_multisampling_early_stop_std: float = 0.04,
    affinity_multisampling_early_stop_patience: int = 1,
    affinity_multisampling_robust_outlier_filter: bool = True,
    affinity_multisampling_robust_outlier_zmax: float = 3.5,
    affinity_multisampling_bootstrap_samples: int = 300,
    max_msa_seqs: int = 8192,
    sampling_steps_affinity: int = 300,
    diffusion_samples_affinity: int = 7,
    cofactor_info: Union[List[Dict], Dict] = None,
    binding_pocket_constraints: Optional[Dict] = None,
    enable_retries: bool = True,
    max_retry_attempts: int = 2,
    retry_delay_base: int = 5,
    subsample_msa: bool = False,
    num_subsampled_msa: int = 1024,
    template_cif_path: Optional[str] = None,
    structure_only: bool = False,
    ptm_modifications: Optional[Dict] = None,
    prediction_timeout_seconds: int = 300,
    enable_msa_cache: bool = True,
    accelerator: str = "gpu",
    devices: int = 1,
    cuda_visible_devices: Optional[str] = None,
    preprocessing_threads: int = 1,
    enable_batch_execution: bool = True,
    use_potentials: bool = False,
    method: Optional[str] = None,
    mutation_steering_config: Optional[Dict[str, Any]] = None,
) -> Tuple[List[Dict], float]:
    """Run screening prediction using Boltz2."""
    results: List[Dict] = []
    protein_drug_filter = st.session_state.get('protein_drug_filter')
    total_prediction_jobs = calculate_filtered_job_count(
        protein_sequences,
        drug_smiles,
        structure_only,
        protein_drug_filter,
    )

    if total_prediction_jobs == 0:
        if protein_drug_filter and protein_drug_filter.get('enabled'):
            st.error("The protein-drug filter excludes all combinations. No predictions will be run.")
        else:
            st.error("No valid protein-drug combinations found.")
        return [], 0.0

    progress_bar = st.progress(0, text="Starting screening prediction...")
    current_job = 0
    existing_results_used = 0
    start_time = time.time()
    wt_sequence_for_msa_cache = next(
        (seq for name, seq in protein_sequences if str(name).strip().upper() == "WT"),
        None
    )
    pending_jobs: List[Dict[str, Any]] = []

    def build_failure_result(
        protein_name: str,
        protein_seq: str,
        drug_name: str,
        drug_smiles_str: str,
        workspace_name: str,
        design_name: str,
        command_params: Dict[str, Any],
        status: str = "Failed - Prediction error",
    ) -> Dict[str, Any]:
        return {
            "protein_name": protein_name,
            "drug_name": "" if structure_only else drug_name,
            "protein_sequence": protein_seq,
            "smiles": "" if structure_only else drug_smiles_str,
            "ic50_um": None,
            "pic50": None,
            "affinity_probability": None,
            "confidence": None,
            "ptm": None,
            "iptm": None,
            "avg_plddt": None,
            "status": status,
            "workspace": workspace_name,
            "design": design_name,
            "cofactor_info": cofactor_info,
            "boltz2_parameters": {
                "use_gpu": use_gpu,
                "accelerator": accelerator,
                "devices": devices,
                "cuda_visible_devices": cuda_visible_devices,
                "preprocessing_threads": preprocessing_threads,
                "recycling_steps": recycling_steps,
                "sampling_steps": sampling_steps,
                "diffusion_samples": diffusion_samples,
                "max_parallel_samples": max_parallel_samples,
                "step_scale": step_scale,
                "affinity_mw_correction": affinity_mw_correction,
                "affinity_consensus_enabled": affinity_consensus_enabled,
                "affinity_consensus_mode": affinity_consensus_mode,
                "affinity_consensus_weight_floor": affinity_consensus_weight_floor,
                "affinity_consensus_entropy_alpha": affinity_consensus_entropy_alpha,
                "external_boltz_patch_enabled": external_boltz_patch_enabled,
                "external_boltz_patch_mode": external_boltz_patch_mode,
                "external_boltz_patch_weight_floor": external_boltz_patch_weight_floor,
                "external_boltz_patch_entropy_alpha": external_boltz_patch_entropy_alpha,
                "external_boltz_patch_uncertainty_penalty": external_boltz_patch_uncertainty_penalty,
                "external_boltz_patch_min_confidence": external_boltz_patch_min_confidence,
                "affinity_multisampling_enabled": affinity_multisampling_enabled,
                "affinity_multisampling_profiles": affinity_multisampling_profiles,
                "affinity_multisampling_settings": affinity_multisampling_settings,
                "affinity_multisampling_refinement_steps": affinity_multisampling_refinement_steps,
                "affinity_multisampling_aggregate_mode": affinity_multisampling_aggregate_mode,
                "affinity_multisampling_apply_aggregate": affinity_multisampling_apply_aggregate,
                "affinity_multisampling_early_stop_enabled": affinity_multisampling_early_stop_enabled,
                "affinity_multisampling_early_stop_min_points": affinity_multisampling_early_stop_min_points,
                "affinity_multisampling_early_stop_delta": affinity_multisampling_early_stop_delta,
                "affinity_multisampling_early_stop_std": affinity_multisampling_early_stop_std,
                "affinity_multisampling_early_stop_patience": affinity_multisampling_early_stop_patience,
                "affinity_multisampling_robust_outlier_filter": affinity_multisampling_robust_outlier_filter,
                "affinity_multisampling_robust_outlier_zmax": affinity_multisampling_robust_outlier_zmax,
                "affinity_multisampling_bootstrap_samples": affinity_multisampling_bootstrap_samples,
                "confidence_target": "balanced",
                "use_potentials": use_potentials,
                "method": method,
                "max_msa_seqs": max_msa_seqs,
                "sampling_steps_affinity": sampling_steps_affinity,
                "diffusion_samples_affinity": diffusion_samples_affinity,
                "subsample_msa": subsample_msa,
                "num_subsampled_msa": num_subsampled_msa,
                "enable_retries": enable_retries,
                "max_retry_attempts": max_retry_attempts,
                "retry_delay_base": retry_delay_base,
                "template_cif_path": template_cif_path,
                "use_cached_msa": command_params.get("use_cached_msa", False),
                "enable_msa_cache": enable_msa_cache,
                "override": command_params["override"],
                "prediction_timeout_seconds": prediction_timeout_seconds,
            },
        }

    ordered_protein_sequences = _order_protein_sequences_for_screening(protein_sequences)

    for protein_name, protein_seq in ordered_protein_sequences:
        combos = [("", "")] if structure_only else drug_smiles
        for drug_name, drug_smiles_str in combos:
            if structure_only:
                if not should_evaluate_protein_drug_pair(protein_name, None, protein_drug_filter):
                    continue
            else:
                if not should_evaluate_protein_drug_pair(protein_name, drug_name, protein_drug_filter):
                    continue

            current_job += 1
            elapsed_time = time.time() - start_time
            if current_job > 1:
                avg_time_per_job = elapsed_time / (current_job - 1)
                remaining_jobs = max(total_prediction_jobs - current_job, 0)
                eta_seconds = remaining_jobs * avg_time_per_job
                eta_text = f"ETA: {eta_seconds / 60:.1f} min"
            else:
                eta_text = "ETA: Calculating..."

            progress = current_job / total_prediction_jobs
            job_label = protein_name if structure_only else f"{protein_name} + {drug_name}"
            progress_bar.progress(progress, text=f"Processing {job_label} ({current_job}/{total_prediction_jobs}) - {eta_text}")

            command_params = {
                "use_gpu": use_gpu,
                "accelerator": accelerator,
                "devices": devices,
                "cuda_visible_devices": cuda_visible_devices,
                "preprocessing_threads": preprocessing_threads,
                "override": not use_existing_results,
                "recycling_steps": recycling_steps,
                "sampling_steps": sampling_steps,
                "diffusion_samples": diffusion_samples,
                "max_parallel_samples": max_parallel_samples,
                "step_scale": step_scale,
                "affinity_mw_correction": affinity_mw_correction,
                "affinity_consensus_enabled": affinity_consensus_enabled,
                "affinity_consensus_mode": affinity_consensus_mode,
                "affinity_consensus_weight_floor": affinity_consensus_weight_floor,
                "affinity_consensus_entropy_alpha": affinity_consensus_entropy_alpha,
                "external_boltz_patch_enabled": external_boltz_patch_enabled,
                "external_boltz_patch_mode": external_boltz_patch_mode,
                "external_boltz_patch_weight_floor": external_boltz_patch_weight_floor,
                "external_boltz_patch_entropy_alpha": external_boltz_patch_entropy_alpha,
                "external_boltz_patch_uncertainty_penalty": external_boltz_patch_uncertainty_penalty,
                "external_boltz_patch_min_confidence": external_boltz_patch_min_confidence,
                "affinity_multisampling_enabled": affinity_multisampling_enabled,
                "affinity_multisampling_profiles": affinity_multisampling_profiles,
                "affinity_multisampling_settings": affinity_multisampling_settings,
                "affinity_multisampling_refinement_steps": affinity_multisampling_refinement_steps,
                "affinity_multisampling_aggregate_mode": affinity_multisampling_aggregate_mode,
                "affinity_multisampling_apply_aggregate": affinity_multisampling_apply_aggregate,
                "affinity_multisampling_early_stop_enabled": affinity_multisampling_early_stop_enabled,
                "affinity_multisampling_early_stop_min_points": affinity_multisampling_early_stop_min_points,
                "affinity_multisampling_early_stop_delta": affinity_multisampling_early_stop_delta,
                "affinity_multisampling_early_stop_std": affinity_multisampling_early_stop_std,
                "affinity_multisampling_early_stop_patience": affinity_multisampling_early_stop_patience,
                "affinity_multisampling_robust_outlier_filter": affinity_multisampling_robust_outlier_filter,
                "affinity_multisampling_robust_outlier_zmax": affinity_multisampling_robust_outlier_zmax,
                "affinity_multisampling_bootstrap_samples": affinity_multisampling_bootstrap_samples,
                "use_potentials": use_potentials,
                "method": method,
                "max_msa_seqs": max_msa_seqs,
                "sampling_steps_affinity": sampling_steps_affinity,
                "diffusion_samples_affinity": diffusion_samples_affinity,
                "enable_retries": enable_retries,
                "max_retry_attempts": max_retry_attempts,
                "retry_delay_base": retry_delay_base,
                "subsample_msa": subsample_msa,
                "num_subsampled_msa": num_subsampled_msa,
                "template_cif_path": template_cif_path,
                "mutation_steering_config": mutation_steering_config,
                "binding_pocket_constraints": _resolve_mutation_steering_constraints(
                    protein_name=protein_name,
                    protein_sequence=protein_seq,
                    base_constraints=binding_pocket_constraints,
                    structure_only=structure_only,
                    mutation_steering_config=mutation_steering_config,
                ),
                "cofactor_info": cofactor_info,
                "ptm_modifications": ptm_modifications,
                "structure_only": structure_only,
                "prediction_timeout_seconds": prediction_timeout_seconds,
                "enable_msa_cache": enable_msa_cache,
            }
            if wt_sequence_for_msa_cache:
                command_params["wt_sequence_for_msa_cache"] = wt_sequence_for_msa_cache

            workspace_name, design_name, stable_job_id = _derive_stable_job_identifiers(
                project_name=project_name,
                protein_name=protein_name,
                protein_sequence=protein_seq,
                drug_name="" if structure_only else drug_name,
                smiles="" if structure_only else drug_smiles_str,
                structure_only=structure_only,
                params=command_params,
            )
            boltz_results = None

            if use_existing_results:
                try:
                    deterministic_yaml_name = f"{workspace_name}_{design_name}"
                    existing_yaml_name = _find_existing_results_by_yaml_name(
                        project_name=project_name,
                        yaml_name=deterministic_yaml_name,
                        structure_only=structure_only,
                    )
                    if not existing_yaml_name:
                        existing_yaml_name = find_existing_screening_results(
                            protein_name,
                            drug_name if not structure_only else "",
                            project_name,
                        )
                except Exception:
                    existing_yaml_name = None
                if existing_yaml_name:
                    try:
                        project_dir = os.path.join(RESULTS_DIR, project_name)
                        yaml_filepath = os.path.join(project_dir, f"{existing_yaml_name}.yaml")
                        is_valid, validation_error = validate_boltz_results(yaml_filepath, structure_only=structure_only)
                        if is_valid:
                            boltz_results = get_screening_existing_boltz_results(
                                existing_yaml_name,
                                project_name,
                                binding_pocket_constraints=command_params.get("binding_pocket_constraints"),
                                cofactor_info=command_params.get("cofactor_info"),
                                structure_only=structure_only,
                            )
                            if boltz_results:
                                workspace_name, design_name = existing_yaml_name, existing_yaml_name
                                existing_results_used += 1
                                progress_bar.progress(progress, text=f"Loading existing results for {job_label} ({current_job}/{total_prediction_jobs}) - {eta_text}")
                            else:
                                st.warning(f"Could not parse existing results for {job_label}")
                        else:
                            st.warning(f"Existing results for {job_label} are invalid: {validation_error}")
                    except Exception as exc:
                        st.warning(f"Failed to load existing results for {job_label}: {str(exc)[:100]}")

            if boltz_results:
                result = _create_result_entry(
                    protein_name=protein_name,
                    protein_sequence=protein_seq,
                    drug_name="" if structure_only else drug_name,
                    smiles="" if structure_only else drug_smiles_str,
                    workspace_name=workspace_name,
                    design_name=design_name,
                    boltz_results=boltz_results,
                    params=command_params,
                )
                results.append(result)
            else:
                pending_jobs.append(
                    {
                        "stable_job_id": stable_job_id,
                        "protein_name": protein_name,
                        "protein_sequence": protein_seq,
                        "drug_name": drug_name,
                        "drug_smiles": drug_smiles_str,
                        "workspace_name": workspace_name,
                        "design_name": design_name,
                        "command_params": command_params,
                    }
                )

    new_computations = len(pending_jobs)

    # Seed WT MSA cache once if needed (mutation mode), then batch remaining jobs.
    if pending_jobs:
        has_mutant_jobs = any(
            (wt_sequence_for_msa_cache and job["protein_sequence"] != wt_sequence_for_msa_cache)
            for job in pending_jobs
        )
        if enable_msa_cache and wt_sequence_for_msa_cache and has_mutant_jobs:
            wt_msa_path, wt_cached = resolve_msa_cache_for_prediction(
                protein_sequence=wt_sequence_for_msa_cache,
                wt_sequence=None,
                enable_msa_cache=enable_msa_cache,
            )
            if not wt_cached:
                wt_seed_idx = next(
                    (i for i, job in enumerate(pending_jobs) if job["protein_sequence"] == wt_sequence_for_msa_cache),
                    None,
                )
                if wt_seed_idx is not None:
                    wt_seed_job = pending_jobs.pop(wt_seed_idx)
                    cp = wt_seed_job["command_params"]
                    seed_start = time.time()
                    boltz_results, success, error_message = run_boltz_with_retry(
                        workspace_name=wt_seed_job["workspace_name"],
                        design_name=wt_seed_job["design_name"],
                        protein_sequence=wt_seed_job["protein_sequence"],
                        ligand_smiles="" if structure_only else wt_seed_job["drug_smiles"],
                        project_name=project_name,
                        protein_display_name=wt_seed_job["protein_name"],
                        ligand_display_name="" if structure_only else wt_seed_job["drug_name"],
                        use_gpu=use_gpu,
                        binding_pocket_constraints=cp.get("binding_pocket_constraints"),
                        override=cp["override"],
                        recycling_steps=cp.get("recycling_steps", recycling_steps),
                        sampling_steps=cp.get("sampling_steps", sampling_steps),
                        diffusion_samples=cp.get("diffusion_samples", diffusion_samples),
                        max_parallel_samples=cp.get("max_parallel_samples", max_parallel_samples),
                        step_scale=cp.get("step_scale", step_scale),
                        affinity_mw_correction=cp.get("affinity_mw_correction", affinity_mw_correction),
                        external_boltz_patch_enabled=cp.get("external_boltz_patch_enabled", False),
                        external_boltz_patch_mode=cp.get("external_boltz_patch_mode", "mutation_aware_v2"),
                        external_boltz_patch_weight_floor=cp.get("external_boltz_patch_weight_floor", 0.05),
                        external_boltz_patch_entropy_alpha=cp.get("external_boltz_patch_entropy_alpha", 0.20),
                        external_boltz_patch_uncertainty_penalty=cp.get("external_boltz_patch_uncertainty_penalty", 0.15),
                        external_boltz_patch_min_confidence=cp.get("external_boltz_patch_min_confidence", 0.35),
                        affinity_multisampling_enabled=cp.get("affinity_multisampling_enabled", False),
                        affinity_multisampling_profiles=cp.get("affinity_multisampling_profiles"),
                        affinity_multisampling_settings=cp.get("affinity_multisampling_settings"),
                        affinity_multisampling_refinement_steps=cp.get("affinity_multisampling_refinement_steps"),
                        affinity_multisampling_aggregate_mode=cp.get("affinity_multisampling_aggregate_mode", "median"),
                        affinity_multisampling_apply_aggregate=cp.get("affinity_multisampling_apply_aggregate", True),
                        affinity_multisampling_early_stop_enabled=cp.get("affinity_multisampling_early_stop_enabled", True),
                        affinity_multisampling_early_stop_min_points=cp.get("affinity_multisampling_early_stop_min_points", 2),
                        affinity_multisampling_early_stop_delta=cp.get("affinity_multisampling_early_stop_delta", 0.02),
                        affinity_multisampling_early_stop_std=cp.get("affinity_multisampling_early_stop_std", 0.04),
                        affinity_multisampling_early_stop_patience=cp.get("affinity_multisampling_early_stop_patience", 1),
                        affinity_multisampling_robust_outlier_filter=cp.get("affinity_multisampling_robust_outlier_filter", True),
                        affinity_multisampling_robust_outlier_zmax=cp.get("affinity_multisampling_robust_outlier_zmax", 3.5),
                        affinity_multisampling_bootstrap_samples=cp.get("affinity_multisampling_bootstrap_samples", 300),
                        confidence_target=cp.get("confidence_target", "balanced"),
                        max_msa_seqs=cp.get("max_msa_seqs", max_msa_seqs),
                        sampling_steps_affinity=cp.get("sampling_steps_affinity", sampling_steps_affinity),
                        diffusion_samples_affinity=cp.get("diffusion_samples_affinity", diffusion_samples_affinity),
                        cofactor_info=cofactor_info,
                        enable_retries=cp.get("enable_retries", enable_retries),
                        max_retry_attempts=cp.get("max_retry_attempts", max_retry_attempts),
                        retry_delay_base=cp.get("retry_delay_base", retry_delay_base),
                        subsample_msa=cp.get("subsample_msa", subsample_msa),
                        num_subsampled_msa=cp.get("num_subsampled_msa", num_subsampled_msa),
                        template_cif_path=cp.get("template_cif_path", template_cif_path),
                        structure_only=structure_only,
                        ptm_modifications=ptm_modifications,
                        prediction_timeout_seconds=cp.get("prediction_timeout_seconds", prediction_timeout_seconds),
                        accelerator=cp.get("accelerator", accelerator),
                        devices=cp.get("devices", devices),
                        cuda_visible_devices=cp.get("cuda_visible_devices", cuda_visible_devices),
                        preprocessing_threads=cp.get("preprocessing_threads", preprocessing_threads),
                        use_potentials=cp.get("use_potentials", False),
                        method=cp.get("method"),
                        msa_path=None,
                        use_cached_msa=False,
                        enable_msa_cache=enable_msa_cache,
                    )
                    seed_elapsed = time.time() - seed_start
                    if success and boltz_results:
                        seed_result = _create_result_entry(
                            protein_name=wt_seed_job["protein_name"],
                            protein_sequence=wt_seed_job["protein_sequence"],
                            drug_name="" if structure_only else wt_seed_job["drug_name"],
                            smiles="" if structure_only else wt_seed_job["drug_smiles"],
                            workspace_name=wt_seed_job["workspace_name"],
                            design_name=wt_seed_job["design_name"],
                            boltz_results=boltz_results,
                            params=cp,
                        )
                    else:
                        seed_result = build_failure_result(
                            protein_name=wt_seed_job["protein_name"],
                            protein_seq=wt_seed_job["protein_sequence"],
                            drug_name=wt_seed_job["drug_name"],
                            drug_smiles_str=wt_seed_job["drug_smiles"],
                            workspace_name=wt_seed_job["workspace_name"],
                            design_name=wt_seed_job["design_name"],
                            command_params=cp,
                            status=f"Failed - {str(error_message)[:120]}",
                        )
                    results.append(seed_result)
                    _persist_result_entry(
                        project_name,
                        seed_result,
                        seed_elapsed,
                        binding_pocket_constraints=cp.get("binding_pocket_constraints"),
                        template_cif_path=template_cif_path,
                        boltz_command=_build_boltz_command(
                            f"{wt_seed_job['workspace_name']}_{wt_seed_job['design_name']}.yaml",
                            cp,
                        ),
                    )

        # Resolve per-job MSA usage before batch or individual fallback.
        for job in pending_jobs:
            cp = job["command_params"]
            msa_path, use_cached_msa = resolve_msa_cache_for_prediction(
                protein_sequence=job["protein_sequence"],
                wt_sequence=wt_sequence_for_msa_cache,
                enable_msa_cache=enable_msa_cache,
            )
            cp["use_cached_msa"] = use_cached_msa
            if msa_path:
                cp["msa_path"] = msa_path
                job["msa_path"] = msa_path
            else:
                job["msa_path"] = None

        remaining_jobs = pending_jobs
        fallback_jobs: List[Dict[str, Any]] = []

        if enable_batch_execution and len(remaining_jobs) >= 2:
            if any(bool(job["command_params"].get("affinity_multisampling_enabled", False)) for job in remaining_jobs):
                fallback_jobs = remaining_jobs
                st.info(
                    "Affinity multi-sampling is enabled; running per-job mode to reuse each job's "
                    "cached structure and execute setting sweeps."
                )
            else:
                batch_compat_keys = [
                    "accelerator",
                    "devices",
                    "cuda_visible_devices",
                    "preprocessing_threads",
                    "override",
                    "recycling_steps",
                    "sampling_steps",
                    "diffusion_samples",
                    "max_parallel_samples",
                    "step_scale",
                    "affinity_mw_correction",
                    "external_boltz_patch_enabled",
                    "external_boltz_patch_mode",
                    "external_boltz_patch_weight_floor",
                    "external_boltz_patch_entropy_alpha",
                    "external_boltz_patch_uncertainty_penalty",
                    "external_boltz_patch_min_confidence",
                    "affinity_multisampling_enabled",
                    "affinity_multisampling_profiles",
                    "affinity_multisampling_settings",
                    "affinity_multisampling_refinement_steps",
                    "affinity_multisampling_aggregate_mode",
                    "affinity_multisampling_apply_aggregate",
                    "affinity_multisampling_early_stop_enabled",
                    "affinity_multisampling_early_stop_min_points",
                    "affinity_multisampling_early_stop_delta",
                    "affinity_multisampling_early_stop_std",
                    "affinity_multisampling_early_stop_patience",
                    "affinity_multisampling_robust_outlier_filter",
                    "affinity_multisampling_robust_outlier_zmax",
                    "affinity_multisampling_bootstrap_samples",
                    "use_potentials",
                    "method",
                    "max_msa_seqs",
                    "sampling_steps_affinity",
                    "diffusion_samples_affinity",
                    "subsample_msa",
                    "num_subsampled_msa",
                ]
                unique_param_signatures = {
                    tuple(_normalize_for_hash(job["command_params"].get(k)) for k in batch_compat_keys)
                    for job in remaining_jobs
                }
                if len(unique_param_signatures) > 1:
                    fallback_jobs = remaining_jobs
                    st.info(
                        "Adaptive sampling produced heterogeneous per-job settings; "
                        "running per-job mode instead of single batch invocation."
                    )
                else:
                    project_dir = os.path.join(RESULTS_DIR, project_name)
                    batch_id = f"batch_inputs_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
                    batch_input_dir = os.path.join(project_dir, batch_id)
                    os.makedirs(batch_input_dir, exist_ok=True)

                    any_job_needs_msa_server = any(not j["command_params"].get("use_cached_msa", False) for j in remaining_jobs)
                    for job in remaining_jobs:
                        cp = job["command_params"]
                        create_screening_boltz_yaml(
                            workspace_name=job["workspace_name"],
                            design_name=job["design_name"],
                            protein_sequence=job["protein_sequence"],
                            ligand_smiles="" if structure_only else job["drug_smiles"],
                            project_name=project_name,
                            binding_pocket_constraints=cp.get("binding_pocket_constraints"),
                            cofactor_info=cofactor_info,
                            template_cif_path=template_cif_path,
                            structure_only=structure_only,
                            ptm_modifications=ptm_modifications,
                            msa_path=job.get("msa_path"),
                            target_directory=batch_input_dir,
                        )

                    batch_start = time.time()
                    batch_command_params = remaining_jobs[0]["command_params"] if remaining_jobs else {}
                    try:
                        utils.run_boltz_batch_prediction(
                            input_path=batch_input_dir,
                            working_dir=project_dir,
                            use_msa_server=any_job_needs_msa_server,
                            accelerator=batch_command_params.get("accelerator", accelerator),
                            devices=batch_command_params.get("devices", devices),
                            cuda_visible_devices=batch_command_params.get("cuda_visible_devices", cuda_visible_devices),
                            preprocessing_threads=batch_command_params.get("preprocessing_threads", preprocessing_threads),
                            override=batch_command_params.get("override", not use_existing_results),
                            recycling_steps=batch_command_params.get("recycling_steps", recycling_steps),
                            sampling_steps=batch_command_params.get("sampling_steps", sampling_steps),
                            diffusion_samples=batch_command_params.get("diffusion_samples", diffusion_samples),
                            max_parallel_samples=batch_command_params.get("max_parallel_samples", max_parallel_samples),
                            step_scale=batch_command_params.get("step_scale", step_scale),
                            affinity_mw_correction=batch_command_params.get("affinity_mw_correction", affinity_mw_correction),
                            use_potentials=batch_command_params.get("use_potentials", False),
                            method=batch_command_params.get("method"),
                            external_boltz_patch_enabled=batch_command_params.get("external_boltz_patch_enabled", False),
                            external_boltz_patch_mode=batch_command_params.get("external_boltz_patch_mode", "mutation_aware_v2"),
                            external_boltz_patch_weight_floor=batch_command_params.get("external_boltz_patch_weight_floor", 0.05),
                            external_boltz_patch_entropy_alpha=batch_command_params.get("external_boltz_patch_entropy_alpha", 0.20),
                            external_boltz_patch_uncertainty_penalty=batch_command_params.get("external_boltz_patch_uncertainty_penalty", 0.15),
                            external_boltz_patch_min_confidence=batch_command_params.get("external_boltz_patch_min_confidence", 0.35),
                            max_msa_seqs=batch_command_params.get("max_msa_seqs", max_msa_seqs),
                            sampling_steps_affinity=batch_command_params.get("sampling_steps_affinity", sampling_steps_affinity),
                            diffusion_samples_affinity=batch_command_params.get("diffusion_samples_affinity", diffusion_samples_affinity),
                            subsample_msa=batch_command_params.get("subsample_msa", subsample_msa),
                            num_subsampled_msa=batch_command_params.get("num_subsampled_msa", num_subsampled_msa),
                            timeout=max(
                                batch_command_params.get("prediction_timeout_seconds", prediction_timeout_seconds) * max(1, len(remaining_jobs)),
                                batch_command_params.get("prediction_timeout_seconds", prediction_timeout_seconds),
                            ),
                        )
                        batch_elapsed = time.time() - batch_start
                        batch_output_root = os.path.join(project_dir, f"boltz_results_{batch_id}")
                        batch_cmd = _build_boltz_batch_command(batch_id, batch_command_params, any_job_needs_msa_server)

                        for job in remaining_jobs:
                            yaml_name = f"{job['workspace_name']}_{job['design_name']}"
                            parsed = utils.parse_boltz_results_from_output_root(
                                batch_output_root,
                                yaml_name,
                                structure_only=structure_only,
                            )
                            if parsed:
                                result = _create_result_entry(
                                    protein_name=job["protein_name"],
                                    protein_sequence=job["protein_sequence"],
                                    drug_name="" if structure_only else job["drug_name"],
                                    smiles="" if structure_only else job["drug_smiles"],
                                    workspace_name=job["workspace_name"],
                                    design_name=job["design_name"],
                                    boltz_results=parsed,
                                    params=job["command_params"],
                                )
                                results.append(result)
                                _persist_result_entry(
                                    project_name,
                                    result,
                                    batch_elapsed / max(len(remaining_jobs), 1),
                                    binding_pocket_constraints=job["command_params"].get("binding_pocket_constraints"),
                                    template_cif_path=template_cif_path,
                                    boltz_command=batch_cmd,
                                )
                            else:
                                fallback_jobs.append(job)
                    except Exception as exc:
                        st.warning(f"Batch execution failed, falling back to per-job execution: {str(exc)[:180]}")
                        fallback_jobs = remaining_jobs
        else:
            fallback_jobs = remaining_jobs

        # Fallback: run any unresolved jobs one-by-one with retry.
        for job in fallback_jobs:
            cp = job["command_params"]
            job_start_time = time.time()
            boltz_results, success, error_message = run_boltz_with_retry(
                workspace_name=job["workspace_name"],
                design_name=job["design_name"],
                protein_sequence=job["protein_sequence"],
                ligand_smiles="" if structure_only else job["drug_smiles"],
                project_name=project_name,
                protein_display_name=job["protein_name"],
                ligand_display_name="" if structure_only else job["drug_name"],
                use_gpu=use_gpu,
                binding_pocket_constraints=cp.get("binding_pocket_constraints"),
                override=cp["override"],
                recycling_steps=cp.get("recycling_steps", recycling_steps),
                sampling_steps=cp.get("sampling_steps", sampling_steps),
                diffusion_samples=cp.get("diffusion_samples", diffusion_samples),
                max_parallel_samples=cp.get("max_parallel_samples", max_parallel_samples),
                step_scale=cp.get("step_scale", step_scale),
                affinity_mw_correction=cp.get("affinity_mw_correction", affinity_mw_correction),
                external_boltz_patch_enabled=cp.get("external_boltz_patch_enabled", False),
                external_boltz_patch_mode=cp.get("external_boltz_patch_mode", "mutation_aware_v2"),
                external_boltz_patch_weight_floor=cp.get("external_boltz_patch_weight_floor", 0.05),
                external_boltz_patch_entropy_alpha=cp.get("external_boltz_patch_entropy_alpha", 0.20),
                external_boltz_patch_uncertainty_penalty=cp.get("external_boltz_patch_uncertainty_penalty", 0.15),
                external_boltz_patch_min_confidence=cp.get("external_boltz_patch_min_confidence", 0.35),
                affinity_multisampling_enabled=cp.get("affinity_multisampling_enabled", False),
                affinity_multisampling_profiles=cp.get("affinity_multisampling_profiles"),
                affinity_multisampling_settings=cp.get("affinity_multisampling_settings"),
                affinity_multisampling_refinement_steps=cp.get("affinity_multisampling_refinement_steps"),
                affinity_multisampling_aggregate_mode=cp.get("affinity_multisampling_aggregate_mode", "median"),
                affinity_multisampling_apply_aggregate=cp.get("affinity_multisampling_apply_aggregate", True),
                affinity_multisampling_early_stop_enabled=cp.get("affinity_multisampling_early_stop_enabled", True),
                affinity_multisampling_early_stop_min_points=cp.get("affinity_multisampling_early_stop_min_points", 2),
                affinity_multisampling_early_stop_delta=cp.get("affinity_multisampling_early_stop_delta", 0.02),
                affinity_multisampling_early_stop_std=cp.get("affinity_multisampling_early_stop_std", 0.04),
                affinity_multisampling_early_stop_patience=cp.get("affinity_multisampling_early_stop_patience", 1),
                affinity_multisampling_robust_outlier_filter=cp.get("affinity_multisampling_robust_outlier_filter", True),
                affinity_multisampling_robust_outlier_zmax=cp.get("affinity_multisampling_robust_outlier_zmax", 3.5),
                affinity_multisampling_bootstrap_samples=cp.get("affinity_multisampling_bootstrap_samples", 300),
                confidence_target=cp.get("confidence_target", "balanced"),
                max_msa_seqs=cp.get("max_msa_seqs", max_msa_seqs),
                sampling_steps_affinity=cp.get("sampling_steps_affinity", sampling_steps_affinity),
                diffusion_samples_affinity=cp.get("diffusion_samples_affinity", diffusion_samples_affinity),
                cofactor_info=cofactor_info,
                enable_retries=cp.get("enable_retries", enable_retries),
                max_retry_attempts=cp.get("max_retry_attempts", max_retry_attempts),
                retry_delay_base=cp.get("retry_delay_base", retry_delay_base),
                subsample_msa=cp.get("subsample_msa", subsample_msa),
                num_subsampled_msa=cp.get("num_subsampled_msa", num_subsampled_msa),
                template_cif_path=cp.get("template_cif_path", template_cif_path),
                structure_only=structure_only,
                ptm_modifications=ptm_modifications,
                prediction_timeout_seconds=cp.get("prediction_timeout_seconds", prediction_timeout_seconds),
                accelerator=cp.get("accelerator", accelerator),
                devices=cp.get("devices", devices),
                cuda_visible_devices=cp.get("cuda_visible_devices", cuda_visible_devices),
                preprocessing_threads=cp.get("preprocessing_threads", preprocessing_threads),
                use_potentials=cp.get("use_potentials", False),
                method=cp.get("method"),
                msa_path=job.get("msa_path"),
                use_cached_msa=cp.get("use_cached_msa", False),
                enable_msa_cache=cp.get("enable_msa_cache", enable_msa_cache),
            )
            job_elapsed = time.time() - job_start_time
            if success and boltz_results:
                result = _create_result_entry(
                    protein_name=job["protein_name"],
                    protein_sequence=job["protein_sequence"],
                    drug_name="" if structure_only else job["drug_name"],
                    smiles="" if structure_only else job["drug_smiles"],
                    workspace_name=job["workspace_name"],
                    design_name=job["design_name"],
                    boltz_results=boltz_results,
                    params=cp,
                )
            else:
                result = build_failure_result(
                    protein_name=job["protein_name"],
                    protein_seq=job["protein_sequence"],
                    drug_name=job["drug_name"],
                    drug_smiles_str=job["drug_smiles"],
                    workspace_name=job["workspace_name"],
                    design_name=job["design_name"],
                    command_params=cp,
                    status=f"Failed - {str(error_message)[:120]}",
                )
            results.append(result)
            _persist_result_entry(
                project_name,
                result,
                job_elapsed,
                binding_pocket_constraints=cp.get("binding_pocket_constraints"),
                template_cif_path=template_cif_path,
                boltz_command=_build_boltz_command(
                    f"{job['workspace_name']}_{job['design_name']}.yaml",
                    cp,
                ),
            )

    computation_time = time.time() - start_time
    progress_bar.empty()

    if use_existing_results and (existing_results_used > 0 or new_computations > 0):
        retry_info = ""
        if enable_retries:
            retry_info = f" (with {max_retry_attempts} retry attempts)"
        st.info(f"Screening summary: {existing_results_used} existing results loaded, {new_computations} new computations performed{retry_info}")
    elif new_computations > 0:
        retry_info = ""
        if enable_retries:
            retry_info = f" (with {max_retry_attempts} retry attempts)"
        st.info(f"Screening summary: {new_computations} new computations performed{retry_info}")

    return results, computation_time

def display_results_table(results: List[Dict]):
    """
    Display results in an interactive table with sorting and filtering.
    
    Args:
        results: List of prediction result dictionaries
    """
    if not results:
        st.info("No results to display.")
        return
    
    # Deduplicate results before processing
    original_count = len(results)
    deduplicated_results = deduplicate_results(results)
    deduplicated_count = len(deduplicated_results)
    
    # Show deduplication info if duplicates were found
    if deduplicated_count < original_count:
        st.warning(f"Found {original_count - deduplicated_count} duplicate entries. Kept {deduplicated_count} unique entries (most complete and recent).")
    
    # Convert to DataFrame
    df = pd.DataFrame(deduplicated_results)
    
    # Add color coding for status
    def color_status(val):
        if val == "Success":
            return "background-color: lightgreen"
        elif val == "Failed":
            return "background-color: lightcoral"
        else:
            return "background-color: lightyellow"

    def _to_float_local(value):
        try:
            if value is None:
                return None
            parsed = float(value)
            if np.isfinite(parsed):
                return parsed
            return None
        except Exception:
            return None

    has_consensus_values = (
        "affinity_pred_value_consensus" in df.columns
        and pd.to_numeric(df["affinity_pred_value_consensus"], errors="coerce").notna().any()
    )
    has_raw_values = (
        "affinity_pred_value_raw" in df.columns
        and pd.to_numeric(df["affinity_pred_value_raw"], errors="coerce").notna().any()
    )
    default_affinity_source_index = 1 if has_consensus_values else 0
    affinity_value_source = st.radio(
        "Affinity Value Source for Tables",
        options=[
            "Raw single-setting output",
            "Consensus output (multi-sampling aggregate when available)",
        ],
        index=default_affinity_source_index,
        horizontal=True,
        help=(
            "Controls which affinity estimate is shown in IC50/pIC50 tables. "
            "Consensus uses the configured aggregate (for example full median) when present. "
            "Rows without that source fall back to available values."
        ),
    )
    source_mode = "consensus" if affinity_value_source.startswith("Consensus") else "raw"

    display_ic50_values = []
    display_pic50_values = []
    display_prob_values = []
    display_source_labels = []
    for row in df.to_dict(orient="records"):
        raw_pred = _to_float_local(row.get("affinity_pred_value_raw"))
        consensus_pred = _to_float_local(row.get("affinity_pred_value_consensus"))
        raw_prob = _to_float_local(row.get("affinity_probability_raw"))
        consensus_prob = _to_float_local(row.get("affinity_probability_consensus"))
        fallback_pic50 = _to_float_local(row.get("pic50"))
        fallback_ic50 = _to_float_local(row.get("ic50_um"))
        fallback_prob = _to_float_local(row.get("affinity_probability"))

        selected_pred = None
        selected_prob = None
        selected_label = "fallback"
        if source_mode == "consensus" and consensus_pred is not None:
            selected_pred = consensus_pred
            selected_prob = consensus_prob
            selected_label = "consensus"
        elif source_mode == "raw" and raw_pred is not None:
            selected_pred = raw_pred
            selected_prob = raw_prob
            selected_label = "raw"
        elif consensus_pred is not None:
            selected_pred = consensus_pred
            selected_prob = consensus_prob
            selected_label = "consensus (fallback)"
        elif raw_pred is not None:
            selected_pred = raw_pred
            selected_prob = raw_prob
            selected_label = "raw (fallback)"

        if selected_pred is not None:
            row_ic50 = 10 ** float(selected_pred)
            row_pic50 = 6.0 - float(selected_pred)
        else:
            row_pic50 = fallback_pic50
            row_ic50 = fallback_ic50

        if selected_prob is None:
            selected_prob = fallback_prob

        display_ic50_values.append(row_ic50)
        display_pic50_values.append(row_pic50)
        display_prob_values.append(selected_prob)
        display_source_labels.append(selected_label)

    df.loc[:, "ic50_um"] = display_ic50_values
    df.loc[:, "pic50"] = display_pic50_values
    df.loc[:, "affinity_probability"] = display_prob_values
    df.loc[:, "affinity_value_source"] = display_source_labels
    
    # Create summary table with IC50 for each drug and protein combination
    if len(deduplicated_results) > 0:
        # Display summary table and screening prediction summary side by side
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.subheader(":material/analytics: Drug Screening Summary")
            
            # Display summary statistics in 2 columns
            col2a, col2b = st.columns(2)
            
            with col2a:
                total_jobs = len(deduplicated_results)
                successful = len([r for r in deduplicated_results if r["status"] == "Success"])
                st.metric("Total Prediction Jobs", total_jobs)
                st.metric("Successful", successful)
            
            with col2b:
                failed = len([r for r in deduplicated_results if r["status"] != "Success"])
                st.metric("Failed", failed)
                
                if successful > 0:
                    avg_pic50 = np.mean([r["pic50"] for r in deduplicated_results if r["pic50"] is not None])
                    st.metric("Avg pIC50", f"{avg_pic50:.2f}")
            
            # Add computation time if available in session state
            if hasattr(st.session_state, 'last_computation_time') and st.session_state.last_computation_time:
                computation_time = st.session_state.last_computation_time
                hours = int(computation_time // 3600)
                minutes = int((computation_time % 3600) // 60)
                seconds = int(computation_time % 60)
                
                if hours > 0:
                    time_str = f"{hours}h {minutes}m {seconds}s"
                elif minutes > 0:
                    time_str = f"{minutes}m {seconds}s"
                else:
                    time_str = f"{seconds}s"
                
                st.metric("Computation Time", time_str)
            elif hasattr(st.session_state, 'loaded_project_data') and st.session_state.loaded_project_data and st.session_state.loaded_project_data.get('computation_time_seconds'):
                # Try to get computation time from loaded project data
                computation_time = st.session_state.loaded_project_data.get('computation_time_seconds')
                if computation_time:
                    hours = int(computation_time // 3600)
                    minutes = int((computation_time % 3600) // 60)
                    seconds = int(computation_time % 60)
                    
                    if hours > 0:
                        time_str = f"{hours}h {minutes}m {seconds}s"
                    elif minutes > 0:
                        time_str = f"{minutes}m {seconds}s"
                    else:
                        time_str = f"{seconds}s"
                    
                    st.metric("Computation Time", time_str)
                else:
                    st.metric("Computation Time", "N/A")
            else:
                # Show N/A when computation time is not available (backward compatibility)
                st.metric("Computation Time", "N/A")
        
        with col2:
            st.subheader(":material/table_chart: IC50 Summary Table (μM)")
            
            # Create pivot table for IC50 values
            summary_df = df.pivot(index="protein_name", columns="drug_name", values="ic50_um")
            # Round values to 4 decimal places for summary table
            summary_df = summary_df.round(4)
            # Gradient highlight: green (min) to red (max) per column
            def highlight_ic50_gradient(s):
                if s.max() != s.min():
                    norm = (s - s.min()) / (s.max() - s.min())
                else:
                    norm = pd.Series([0.5] * len(s), index=s.index)
                # Use RdYlGn palette: green (low IC50, good binding) -> yellow -> red (high IC50, poor binding)
                def color(val, n):
                    if pd.isnull(val):
                        return ''
                    if n <= 0.5:
                        # Green to Yellow (n: 0->0.5) - good to moderate binding
                        r1, g1, b1 = (26, 152, 80)   # green
                        r2, g2, b2 = (254, 224, 139) # yellow
                        t = n * 2  # scale 0-0.5 to 0-1
                        r = int(r1 + (r2 - r1) * t)
                        g = int(g1 + (g2 - g1) * t)
                        b = int(b1 + (b2 - b1) * t)
                    else:
                        # Yellow to Red (n: 0.5->1) - moderate to poor binding
                        r1, g1, b1 = (254, 224, 139) # yellow
                        r2, g2, b2 = (215, 48, 39)   # red
                        t = (n - 0.5) * 2  # scale 0.5-1 to 0-1
                        r = int(r1 + (r2 - r1) * t)
                        g = int(g1 + (g2 - g1) * t)
                        b = int(b1 + (b2 - b1) * t)
                    # Calculate luminance to determine text color
                    # Use relative luminance formula: 0.299*R + 0.587*G + 0.114*B
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = "white" if luminance < 0.5 else "black"
                    return f'background-color: rgb({r},{g},{b}); color: {text_color}'
                return [color(v, n) for v, n in zip(s, norm)]
            styled_summary_df = summary_df.style.apply(highlight_ic50_gradient, axis=0)
            st.dataframe(
                styled_summary_df,
                use_container_width=True,
                hide_index=False,
                column_config={
                    "protein_name": st.column_config.TextColumn("Protein", width="medium"),
                }
            )
    
    # Display results table
    st.subheader(":material/table: Detailed Results")
    
    # Reorder columns to move protein_sequence, SMILES, and status to the end (exclude workspace and design)
    column_order = [
        "protein_name", "drug_name", "ic50_um", "pic50", "affinity_probability", 
        "affinity_value_source",
        "confidence", "mutation_local_consistency_score", "mutation_local_consistency_label",
        "mutation_local_fingerprint_similarity", "mutation_local_disruption_score",
        "ptm", "iptm", "avg_plddt",
        "protein_sequence", "smiles", "status"
    ]
    
    # Filter columns that exist in the dataframe
    existing_columns = [col for col in column_order if col in df.columns]
    df_reordered = df[existing_columns].copy()
    
    # Ensure IC50 is always displayed with 4 decimals in the detailed table
    if "ic50_um" in df_reordered.columns:
        df_reordered.loc[:, "ic50_um"] = df_reordered["ic50_um"].apply(lambda x: round(x, 4) if pd.notnull(x) else x)
    
    # Create styled DataFrame
    styled_df = df_reordered.style.map(color_status, subset=['status'])
    
    # Display with column configuration
    st.dataframe(
        styled_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "protein_name": st.column_config.TextColumn("Protein", width="medium"),
            "drug_name": st.column_config.TextColumn("Drug", width="medium"),
            "ic50_um": st.column_config.NumberColumn("IC50 (μM)", format="%.4f"),
            "pic50": st.column_config.NumberColumn("pIC50", format="%.3f"),
            "affinity_probability": st.column_config.NumberColumn("Affinity Prob", format="%.3f"),
            "affinity_value_source": st.column_config.TextColumn(
                "Affinity Source",
                width="medium",
                help=(
                    "Indicates whether displayed IC50/pIC50 came from raw single-setting output "
                    "or consensus aggregate (for example full median)."
                ),
            ),
            "confidence": st.column_config.NumberColumn("Confidence", format="%.3f"),
            "mutation_local_consistency_score": st.column_config.NumberColumn(
                "Mut Local Consistency",
                format="%.3f",
                help="0-1 score comparing WT vs mutant local interaction consistency for the same drug.",
            ),
            "mutation_local_consistency_label": st.column_config.TextColumn(
                "Consistency Label",
                width="small",
            ),
            "mutation_local_fingerprint_similarity": st.column_config.NumberColumn(
                "FP Similarity",
                format="%.3f",
            ),
            "mutation_local_disruption_score": st.column_config.NumberColumn(
                "Disruption",
                format="%.3f",
                help="Higher values indicate larger WT vs mutant interface disruption.",
            ),
            "ptm": st.column_config.NumberColumn("pTM", format="%.3f"),
            "iptm": st.column_config.NumberColumn("ipTM", format="%.3f"),
            "avg_plddt": st.column_config.NumberColumn("Avg pLDDT", format="%.1f"),
            "protein_sequence": st.column_config.TextColumn("Protein Sequence", width="large", help="Bracketed residues [X] indicate mutations"),
            "smiles": st.column_config.TextColumn("SMILES", width="large"),
            "status": st.column_config.TextColumn("Status", width="medium")
        }
    )
    
    col1, _, col3 = st.columns(3)
    with col1:
        # pIC50 to IC50 Converter
        with st.popover("pIC50 to IC50 (µM) Converter", icon=":material/swap_horiz:"):
            st.markdown("""
            **Convert between pIC50 and IC50 (μM)**
            """)
            col1, col2 = st.columns(2)
            with col1:
                pic50_input = st.text_input("pIC50", key="pic50_to_ic50_input", help="pIC50 = -log10(IC50 [M])", placeholder="e.g. 7.5")
                ic50_result = None
                if pic50_input:
                    try:
                        pic50_val = float(pic50_input)
                        ic50_um = 10 ** (-pic50_val) * 1e6
                        ic50_result = f"→ IC50 = {ic50_um:.4g} μM"
                        st.write(ic50_result)
                    except Exception:
                        st.write("Invalid pIC50 value")
            with col2:
                ic50_input = st.text_input("IC50 (μM)", key="ic50_to_pic50_input", help="IC50 [μM] = 10^(-pIC50) × 1e6", placeholder="e.g. 0.5")
                pic50_result = None
                if ic50_input:
                    try:
                        ic50_val = float(ic50_input)
                        if ic50_val <= 0:
                            raise ValueError
                        pic50 = -math.log10(ic50_val * 1e-6)
                        pic50_result = f"→ pIC50 = {pic50:.4g}"
                        st.write(pic50_result)
                    except Exception:
                        st.write("Invalid IC50 value (must be > 0)")

    with col3:
        # Add rename functionality popup
        with st.popover("Rename Results", icon=":material/edit:", help="Rename protein or drug names and update all associated files"):
            st.write("**Rename protein or drug names in all results**")
            
            # Get unique protein and drug names
            unique_proteins = sorted(df_reordered['protein_name'].unique())
            unique_drugs = sorted(df_reordered['drug_name'].unique())
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Rename Proteins")
                if unique_proteins:
                    old_protein = st.selectbox("Select protein to rename:", unique_proteins, key="rename_protein_select")
                    new_protein = st.text_input("New name:", value=old_protein, key="rename_protein_new")
                    
                    if st.button("Rename Protein", key="rename_protein_btn"):
                        if old_protein != new_protein and new_protein.strip():
                            success = rename_results_in_project(
                                project_name=st.session_state.get('loaded_project_name'),
                                old_name=old_protein,
                                new_name=new_protein.strip(),
                                rename_type="protein",
                                results_dir=RESULTS_DIR
                            )
                            if success:
                                st.success(f"Renamed '{old_protein}' to '{new_protein.strip()}' in all files")
                                # Reload results from disk to reflect changes
                                project_data = load_project_data(st.session_state.get('loaded_project_name'), RESULTS_DIR)
                                if project_data and isinstance(project_data, dict) and 'results' in project_data:
                                    st.session_state.screening_results = deduplicate_results(project_data['results'])
                                st.rerun()
                            else:
                                st.error("Failed to rename protein. Check console for details.")
                        else:
                            st.warning("Please enter a different name")
                else:
                    st.info("No proteins found")
            
            with col2:
                st.subheader("Rename Drugs")
                if unique_drugs:
                    old_drug = st.selectbox("Select drug to rename:", unique_drugs, key="rename_drug_select")
                    new_drug = st.text_input("New name:", value=old_drug, key="rename_drug_new")
                    
                    if st.button("Rename Drug", key="rename_drug_btn"):
                        if old_drug != new_drug and new_drug.strip():
                            success = rename_results_in_project(
                                project_name=st.session_state.get('loaded_project_name'),
                                old_name=old_drug,
                                new_name=new_drug.strip(),
                                rename_type="drug",
                                results_dir=RESULTS_DIR
                            )
                            if success:
                                st.success(f"Renamed '{old_drug}' to '{new_drug.strip()}' in all files")
                                # Reload results from disk to reflect changes
                                project_data = load_project_data(st.session_state.get('loaded_project_name'), RESULTS_DIR)
                                if project_data and isinstance(project_data, dict) and 'results' in project_data:
                                    st.session_state.screening_results = deduplicate_results(project_data['results'])
                                st.rerun()
                            else:
                                st.error("Failed to rename drug. Check console for details.")
                        else:
                            st.warning("Please enter a different name")
                else:
                    st.info("No drugs found")
            
            st.write("**Note:** This will update all screening_results_*.json files and the project_metadata.json file in the project folder.")

def force_update_project_from_boltz_results(project_name: str, results_dir: str):
    """
    Force update project_metadata.json by scanning all boltz_results_screening_* folders
    for valid results and extracting data from affinity_*.json and *.pdb files.
    
    Args:
        project_name (str): Name of the project to update
        results_dir (str): Path to the directory containing project folders
    """
    try:
        project_dir = os.path.join(results_dir, project_name)
        if not os.path.exists(project_dir):
            st.error(f"Project directory not found: {project_dir}")
            return
        
        # Find all boltz_results_screening_* directories
        boltz_dirs = []
        for item in os.listdir(project_dir):
            if item.startswith("boltz_results_screening_") and os.path.isdir(os.path.join(project_dir, item)):
                boltz_dirs.append(item)
        
        if not boltz_dirs:
            st.warning(f"No boltz_results_screening_* directories found in project: {project_name}")
            return
        
        all_results = []
        processed_count = 0
        
        for boltz_dir in boltz_dirs:
            boltz_path = os.path.join(project_dir, boltz_dir)
            predictions_path = os.path.join(boltz_path, "predictions")
            
            if not os.path.exists(predictions_path):
                continue
            
            # Look for subdirectories in predictions/
            for pred_subdir in os.listdir(predictions_path):
                pred_subdir_path = os.path.join(predictions_path, pred_subdir)
                if not os.path.isdir(pred_subdir_path):
                    continue
                
                # Look for affinity_*.json files
                affinity_files = glob.glob(os.path.join(pred_subdir_path, "affinity_*.json"))
                
                for affinity_file in affinity_files:
                    # Check if there are also PDB files in the same directory
                    pdb_files = glob.glob(os.path.join(pred_subdir_path, "*.pdb"))
                    
                    if pdb_files:  # Valid results found
                        try:
                            # Parse the affinity results
                            with open(affinity_file, 'r') as f:
                                affinity_data = json.load(f)
                            
                            # Extract design information from directory structure
                            # boltz_results_screening_{timestamp}_{protein}_{drug}
                            dir_parts = boltz_dir.replace("boltz_results_screening_", "").split("_")
                            if len(dir_parts) >= 3:
                                timestamp = dir_parts[0]
                                protein_name = "_".join(dir_parts[1:-1])
                                drug_name = dir_parts[-1]
                                
                                # Create result entry similar to normal screening processing
                                result_entry = {
                                    "protein": protein_name,
                                    "drug": drug_name,
                                    "status": "Success",
                                    "timestamp": timestamp,
                                    "workspace": f"screening_{timestamp}",
                                    "design": pred_subdir,
                                    "boltz_dir": boltz_dir,
                                    "pdb_files": [os.path.basename(pdb) for pdb in pdb_files],
                                }
                                
                                # Add affinity data if available
                                if isinstance(affinity_data, dict):
                                    if "boltz_affinity" in affinity_data:
                                        result_entry["boltz_affinity"] = affinity_data["boltz_affinity"]
                                    if "boltz_affinity_confidence" in affinity_data:
                                        result_entry["boltz_affinity_confidence"] = affinity_data["boltz_affinity_confidence"]
                                
                                # Look for confidence file as well
                                confidence_files = glob.glob(os.path.join(pred_subdir_path, "confidence_*.json"))
                                if confidence_files:
                                    try:
                                        with open(confidence_files[0], 'r') as f:
                                            confidence_data = json.load(f)
                                            if isinstance(confidence_data, dict) and "confidences" in confidence_data:
                                                result_entry["confidences"] = confidence_data["confidences"]
                                    except Exception:
                                        pass  # Skip confidence data if it can't be read
                                
                                all_results.append(result_entry)
                                processed_count += 1
                        
                        except Exception as e:
                            st.warning(f"Error processing {affinity_file}: {str(e)}")
                            continue
        
        if all_results:
            # Save as a new screening results file
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"screening_results_force_update_{timestamp}.json"
            filepath = os.path.join(project_dir, filename)
            
            results_with_metadata = {
                "computation_time_seconds": None,
                "timestamp": timestamp,
                "force_update": True,
                "source": "boltz_results_screening_directories",
                "results": all_results
            }
            
            with open(filepath, 'w') as f:
                json.dump(results_with_metadata, f, indent=4, default=str)
            
            st.success(f"Force update completed: Found and processed {processed_count} valid results from {len(boltz_dirs)} boltz directories")
            st.info(f"Results saved to: {filename}")
        else:
            st.warning(f"No valid results found in boltz_results_screening_* directories for project: {project_name}")
    
    except Exception as e:
        st.error(f"Error during force update: {str(e)}")

def concatenate_screening_results_to_metadata():
    """
    Concatenates all screening_results_*.json files to update project_metadata.json files.
    If there are conflicting entries, keeps the latest by timestamp.
    """
    try:
        # Find all screening results files
        screening_files = glob.glob("**/screening_results_*.json", recursive=True)
        
        if not screening_files:
            return
        
        # Group screening files by project directory
        project_screenings = {}
        for screening_file in screening_files:
            project_dir = os.path.dirname(screening_file)
            if project_dir not in project_screenings:
                project_screenings[project_dir] = []
            project_screenings[project_dir].append(screening_file)
        
        # Process each project directory
        for project_dir, screening_files_list in project_screenings.items():
            metadata_file = os.path.join(project_dir, "project_metadata.json")
            
            # Load existing metadata or create new structure
            if os.path.exists(metadata_file):
                try:
                    with open(metadata_file, 'r') as f:
                        metadata = json.load(f)
                    
                    # Handle old format where metadata is a list instead of dict
                    if isinstance(metadata, list):
                        old_results = metadata
                        metadata = {
                            "project_name": os.path.basename(project_dir),
                            "created_date": datetime.now().isoformat(),
                            "last_updated": datetime.now().isoformat(),
                            "total_results": len(old_results),
                            "successful_results": len([r for r in old_results if isinstance(r, dict) and r.get("status") == "Success"]),
                            "failed_results": len([r for r in old_results if isinstance(r, dict) and r.get("status") != "Success"]),
                            "computation_time_seconds": 0.0,
                            "results": old_results
                        }
                except (json.JSONDecodeError, FileNotFoundError):
                    metadata = {
                        "project_name": os.path.basename(project_dir),
                        "created_date": datetime.now().isoformat(),
                        "last_updated": datetime.now().isoformat(),
                        "total_results": 0,
                        "successful_results": 0,
                        "failed_results": 0,
                        "computation_time_seconds": 0.0,
                        "results": []
                    }
            else:
                metadata = {
                    "project_name": os.path.basename(project_dir),
                    "created_date": datetime.now().isoformat(),
                    "last_updated": datetime.now().isoformat(),
                    "total_results": 0,
                    "successful_results": 0,
                    "failed_results": 0,
                    "computation_time_seconds": 0.0,
                    "results": []
                }
            
            # Collect all results with timestamps
            all_results = {}  # Key: unique identifier, Value: (result, timestamp)
            total_computation_time = metadata.get("computation_time_seconds", 0.0) or 0.0
            latest_timestamp = None
            
            # First, add existing results from metadata with their timestamps
            existing_results = metadata.get("results", [])
            for result in existing_results:
                # Skip if result is not a dictionary
                if not isinstance(result, dict):
                    continue
                    
                # Create unique identifier based on protein_name, drug_name, and design
                protein_name = result.get("protein_name", "")
                drug_name = result.get("drug_name", "")
                design = result.get("design", "")
                
                unique_key = f"{protein_name}_{drug_name}_{design}"
                
                # Try to get timestamp from result, or use a very old timestamp as fallback
                result_timestamp_str = result.get("timestamp", "")
                try:
                    if result_timestamp_str:
                        result_datetime = datetime.fromisoformat(result_timestamp_str.replace('Z', '+00:00'))
                    else:
                        # Use a very old timestamp for existing results without timestamps
                        result_datetime = datetime(1900, 1, 1)
                except (ValueError, TypeError):
                    result_datetime = datetime(1900, 1, 1)
                
                all_results[unique_key] = (result, result_datetime)
                
                if latest_timestamp is None or result_datetime > latest_timestamp:
                    latest_timestamp = result_datetime
            
            # Process each screening file
            for screening_file in screening_files_list:
                try:
                    with open(screening_file, 'r') as f:
                        screening_data = json.load(f)
                    
                    # Handle old format where screening_data is a list instead of dict
                    if isinstance(screening_data, list):
                        screening_timestamp = ""
                        screening_results = screening_data
                        screening_computation_time = 0.0
                    else:
                        screening_timestamp = screening_data.get("timestamp", "")
                        screening_results = screening_data.get("results", [])
                        screening_computation_time = screening_data.get("computation_time_seconds", 0.0)
                    
                    # Convert timestamp to datetime for comparison
                    try:
                        screening_datetime = datetime.strptime(screening_timestamp, "%Y%m%d_%H%M%S")
                        if latest_timestamp is None or screening_datetime > latest_timestamp:
                            latest_timestamp = screening_datetime
                    except ValueError:
                        screening_datetime = datetime.now()
                    
                    # Add computation time (handle None values)
                    if screening_computation_time is not None:
                        total_computation_time += screening_computation_time
                    
                    # Process each result in the screening
                    for result in screening_results:
                        # Skip if result is not a dictionary
                        if not isinstance(result, dict):
                            continue
                            
                        # Create unique identifier based on protein_name, drug_name, and design
                        protein_name = result.get("protein_name", "")
                        drug_name = result.get("drug_name", "")
                        design = result.get("design", "")
                        
                        unique_key = f"{protein_name}_{drug_name}_{design}"
                        
                        # Keep the result with the latest timestamp
                        if unique_key not in all_results or screening_datetime > all_results[unique_key][1]:
                            all_results[unique_key] = (result, screening_datetime)
                
                except (json.JSONDecodeError, FileNotFoundError, KeyError) as e:
                    print(f"Error processing screening file {screening_file}: {e}")
                    continue
            
            # Extract results and calculate statistics
            final_results = [result_data[0] for result_data in all_results.values()]
            successful_count = sum(1 for result in final_results if result.get("status", "").startswith("Success"))
            failed_count = len(final_results) - successful_count
            
            # Update metadata
            metadata.update({
                "last_updated": (latest_timestamp or datetime.now()).isoformat(),
                "total_results": len(final_results),
                "successful_results": successful_count,
                "failed_results": failed_count,
                "computation_time_seconds": total_computation_time,
                "results": final_results
            })
            
            # Save updated metadata
            try:
                with open(metadata_file, 'w') as f:
                    json.dump(metadata, f, indent=4, default=str)
            except Exception as e:
                print(f"Error saving metadata file {metadata_file}: {e}")
    
    except Exception as e:
        print(f"Error in concatenate_screening_results_to_metadata: {e}")

def main():
    """Main function for the drug screening page."""
    st.set_page_config(
        page_title="BoltzOmics",
        page_icon=Image.open(os.path.join("static", "boltzomics", "boltzomics_icon.png")),
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # Load CSS
    load_css()
    
    bg_images = glob.glob("static/background_*.jpg")
    if bg_images:
        bg_image = bg_images[0].replace("static/", "app/static/")
    else:
        bg_image = "app/static/background_1.jpg"
    header_html = f'''
    <div class="header-container" style="background-image: url('{bg_image}');">
        <img src="app/static/boltzomics/boltzomics_light.png" style="width: 30vw; height: auto; max-width: 30%; z-index: 3; position: relative; margin-left: -50%; margin-top: 4rem;" alt="Boltzomics Logo">
    </div>
    '''
    st.markdown(header_html, unsafe_allow_html=True)
    
    # Initialize session state for cofactor info
    if 'cofactor_info' not in st.session_state:
        st.session_state.cofactor_info = []
    structure_only = bool(st.session_state.get("structure_only", False))
    
    # Sidebar with logo and navigation
    with st.sidebar:
        # st.image(os.path.join("static", "boltzomics", "boltzomics_light.png"), use_container_width=True)

        st.header(":material/settings: Configuration")
        
        # Computation Settings
        st.subheader("Compute")
        available_gpus = discover_gpu_devices()
        accelerator_choice = st.selectbox(
            "Run On",
            options=["GPU", "CPU"],
            index=0,
            help="Choose GPU for faster runs, or CPU if no GPU is available.",
        )
        accelerator = accelerator_choice.lower()
        use_gpu = accelerator == "gpu"
        gpu_execution_mode = "single"
        queue_worker_count = 1
        queue_gpu_devices: List[str] = []

        cuda_visible_devices: Optional[str] = None
        if use_gpu:
            if len(available_gpus) > 1:
                gpu_execution_mode = st.selectbox(
                    "GPU Usage Mode",
                    options=[
                        "Single GPU",
                        "Multi-GPU Queue (One Job per GPU)",
                    ],
                    index=0,
                    help=(
                        "Single GPU: run one job at a time on one GPU. "
                        "Multi-GPU Queue: run multiple jobs in parallel, one per selected GPU."
                    ),
                )
            gpu_labels = [f"GPU {idx}: {name}" for idx, name in available_gpus]
            if gpu_execution_mode.startswith("Multi-GPU"):
                default_labels = gpu_labels if gpu_labels else []
                selected_gpu_labels = st.multiselect(
                    "GPUs to Use",
                    options=gpu_labels,
                    default=default_labels,
                    help="Pick which GPUs can run jobs in parallel.",
                )
                if not selected_gpu_labels and gpu_labels:
                    selected_gpu_labels = [gpu_labels[0]]
                queue_gpu_devices = [
                    label.split(":", 1)[0].replace("GPU", "").strip()
                    for label in selected_gpu_labels
                ]
                queue_worker_count = max(1, len(queue_gpu_devices))
                # Non-queue fallback path uses the first selected GPU.
                if queue_gpu_devices:
                    cuda_visible_devices = queue_gpu_devices[0]
            else:
                gpu_labels_with_auto = ["Auto (Default CUDA visibility)"] + gpu_labels
                selected_gpu_label = st.selectbox(
                    "GPU Device",
                    options=gpu_labels_with_auto,
                    index=0,
                    help=(
                        "Pick a specific GPU, or choose Auto to let the system decide."
                    ),
                )
                if selected_gpu_label.startswith("GPU "):
                    cuda_visible_devices = selected_gpu_label.split(":", 1)[0].replace("GPU", "").strip()
                    queue_gpu_devices = [cuda_visible_devices]
                elif not available_gpus:
                    st.caption("No local GPU list detected; Boltz will use default CUDA device visibility.")

        boltz_devices = st.number_input(
            "Boltz Devices",
            min_value=1,
            value=1,
            help="How many devices each single Boltz run can use.",
        )
        if gpu_execution_mode.startswith("Multi-GPU"):
            boltz_devices = 1
            st.caption("Multi-GPU queue mode runs one job per GPU.")
        elif use_gpu and cuda_visible_devices not in (None, "", "auto") and boltz_devices > 1:
            st.caption("A specific GPU is selected, so each job uses one device.")
            boltz_devices = 1

        preprocessing_threads = st.number_input(
            "Preprocessing Threads",
            min_value=1,
            value=max(1, min((os.cpu_count() or 1), 8)),
            help="CPU threads used to prepare inputs. Higher can be faster, but uses more CPU.",
        )

        result_reuse_mode = st.selectbox(
            "Existing Results Policy",
            options=[
                "Reuse Valid Results (Recommended)",
                "Recompute All Jobs",
            ],
            index=0,
            help=(
                "Reuse: skip already-computed protein-drug jobs and load saved outputs. "
                "Recompute: ignore saved outputs and run everything again."
            ),
        )
        use_existing_results = result_reuse_mode.startswith("Reuse")
        enable_batch_execution = st.toggle(
            "Batch New Jobs (Single Boltz Invocation)",
            value=True,
            help=(
                "Run new jobs together in one launch to reduce startup overhead. "
                "Automatically falls back to one-by-one if needed."
            ),
        )
        max_parallel_samples = st.number_input(
            "Max Parallel Samples",
            min_value=1,
            value=5,
            help="How many samples to process at once in a run. Higher can be faster but uses more memory.",
        )

        # Sampling Settings
        st.subheader("Structure Quality")
        recycling_steps = 4
        sampling_steps = 300
        diffusion_samples = 1
        step_scale = 1.638
        recycling_steps = st.number_input("Recycling Steps", min_value=1, value=4, help="Extra refinement rounds for structure prediction. More rounds may improve quality but take longer.")
        sampling_steps = st.number_input("Sampling Steps", min_value=1, value=300, help="Number of structure sampling steps. More steps are usually more stable but slower.")
        diffusion_samples = st.number_input("Diffusion Samples", min_value=1, value=1, help="How many structure samples to generate per job. More samples improve robustness but increase runtime.")
        step_scale = st.number_input("Step Scale", value=1.638, format="%.3f", help="Controls exploration vs precision during sampling. Keep default unless you are tuning.")

        # Affinity Prediction Settings
        st.subheader("Affinity Prediction")
        affinity_mw_correction = False
        affinity_multisampling_enabled = False
        affinity_multisampling_profiles: List[Dict[str, int]] = []
        affinity_multisampling_settings: List[str] = []
        affinity_multisampling_refinement_steps: List[int] = []
        affinity_multisampling_aggregate_mode = "median"
        affinity_multisampling_apply_aggregate = True
        affinity_multisampling_early_stop_enabled = True
        affinity_multisampling_early_stop_min_points = 2
        affinity_multisampling_early_stop_delta = 0.02
        affinity_multisampling_early_stop_std = 0.04
        affinity_multisampling_early_stop_patience = 1
        affinity_multisampling_robust_outlier_filter = True
        affinity_multisampling_robust_outlier_zmax = 3.5
        affinity_multisampling_bootstrap_samples = 300
        sampling_steps_affinity = 300
        diffusion_samples_affinity = 7
        affinity_consensus_enabled = False
        affinity_consensus_mode = "weighted"
        affinity_consensus_weight_floor = 0.05
        affinity_consensus_entropy_alpha = 0.20
        external_boltz_patch_enabled = False
        external_boltz_patch_mode = "mutation_aware_v2"
        external_boltz_patch_weight_floor = 0.05
        external_boltz_patch_entropy_alpha = 0.20
        external_boltz_patch_uncertainty_penalty = 0.15
        external_boltz_patch_min_confidence = 0.35

        if structure_only:
            st.caption("Affinity settings are hidden in Structure-only mode.")
        else:
            affinity_mw_correction = st.toggle(
                "Molecular Weight Correction",
                value=False,
                help="Adjusts affinity for molecular size. Leave OFF unless you specifically want this correction.",
            )
            affinity_multisampling_enabled = st.toggle(
                "Enable Affinity Multi-Sampling (Structure Once, Affinity Sweep)",
                value=True,
                help=(
                    "Run structure once, then compute affinity at multiple settings. "
                    "This gives a more robust final score without rerunning structure each time."
                ),
            )

            if affinity_multisampling_enabled:
                st.caption(
                    "Recommended quality mode: run a sweep and use Full Consensus (Median) as the final affinity."
                )
                affinity_multisampling_steps_text = st.text_input(
                    "Multi Sampling Steps (Affinity)",
                    value="200,300,400",
                    help=(
                        "Comma-separated affinity step values. Example: 200,300,400."
                    ),
                )
                affinity_multisampling_steps = _parse_positive_int_list_from_text(
                    affinity_multisampling_steps_text
                )
                if not affinity_multisampling_steps:
                    affinity_multisampling_steps = [200, 300, 400]
                    st.caption("Using default multi-sampling steps: 200, 300, 400.")

                affinity_multisampling_diffusion_text = st.text_input(
                    "Multi Diffusion Samples (Affinity)",
                    value="5,7,9",
                    help=(
                        "Comma-separated diffusion values matching the steps above (same order). "
                        "Use one value (example: 7) to apply to all. Leave empty for automatic values."
                    ),
                )
                multi_diffusion_values: List[int]
                if not str(affinity_multisampling_diffusion_text).strip():
                    multi_diffusion_values = [
                        _auto_diffusion_samples_for_affinity_step(step)
                        for step in affinity_multisampling_steps
                    ]
                    st.caption("Diffusion samples auto-derived from steps.")
                else:
                    parsed_diffusion = _parse_positive_int_list_from_text(affinity_multisampling_diffusion_text)
                    if not parsed_diffusion:
                        multi_diffusion_values = [
                            _auto_diffusion_samples_for_affinity_step(step)
                            for step in affinity_multisampling_steps
                        ]
                        st.caption("Invalid diffusion list; using auto-derived diffusion samples.")
                    elif len(parsed_diffusion) == 1:
                        multi_diffusion_values = [int(parsed_diffusion[0])] * len(affinity_multisampling_steps)
                    else:
                        adjusted_diffusion = [int(v) for v in parsed_diffusion]
                        if len(adjusted_diffusion) < len(affinity_multisampling_steps):
                            adjusted_diffusion = adjusted_diffusion + [adjusted_diffusion[-1]] * (
                                len(affinity_multisampling_steps) - len(adjusted_diffusion)
                            )
                        elif len(adjusted_diffusion) > len(affinity_multisampling_steps):
                            adjusted_diffusion = adjusted_diffusion[: len(affinity_multisampling_steps)]
                        multi_diffusion_values = adjusted_diffusion

                full_consensus_median_enabled = st.toggle(
                    "Use Full Consensus (Median) as Final Affinity (Recommended)",
                    value=True,
                    help=(
                        "If ON, the app uses the median across all executed multi-sampling affinity settings "
                        "as the main reported affinity value. In results tables, you can switch display "
                        "between raw and consensus values via the `Affinity Value Source` control."
                    ),
                )
                if full_consensus_median_enabled:
                    affinity_multisampling_aggregate_mode = "median"
                    affinity_multisampling_apply_aggregate = True
                else:
                    affinity_multisampling_aggregate_mode = "median"
                    affinity_multisampling_apply_aggregate = False

                affinity_multisampling_profiles = [
                    {
                        "sampling_steps_affinity": int(step),
                        "diffusion_samples_affinity": int(diff),
                    }
                    for step, diff in zip(affinity_multisampling_steps, multi_diffusion_values)
                ]
                affinity_multisampling_refinement_steps = [
                    int(profile["sampling_steps_affinity"]) for profile in affinity_multisampling_profiles
                ]
                affinity_multisampling_settings = [
                    f"{int(profile['sampling_steps_affinity'])}x{int(profile['diffusion_samples_affinity'])}"
                    for profile in affinity_multisampling_profiles
                ]
                sampling_steps_affinity = int(affinity_multisampling_profiles[0]["sampling_steps_affinity"])
                diffusion_samples_affinity = int(affinity_multisampling_profiles[0]["diffusion_samples_affinity"])
                st.caption(
                    "Multi-sampling profiles: "
                    + ", ".join(
                        [
                            f"{int(p['sampling_steps_affinity'])}:{int(p['diffusion_samples_affinity'])}"
                            for p in affinity_multisampling_profiles
                        ]
                    )
                )
            else:
                sampling_steps_affinity = st.number_input(
                    "Sampling Steps (Affinity)",
                    min_value=1,
                    value=300,
                    help="Number of sampling steps for single-run affinity prediction.",
                )
                diffusion_samples_affinity = st.number_input(
                    "Diffusion Samples (Affinity)",
                    min_value=1,
                    value=7,
                    help="Number of diffusion samples for single-run affinity prediction.",
                )

            # Keep the default user-facing affinity workflow focused on
            # multi-sampling with consensus aggregation.
            affinity_consensus_enabled = False
            external_boltz_patch_enabled = False

            if not affinity_multisampling_enabled:
                affinity_multisampling_profiles = []
                affinity_multisampling_settings = []
                affinity_multisampling_refinement_steps = []
        if affinity_multisampling_enabled and enable_batch_execution:
            st.caption(
                "Batch mode will be bypassed when affinity multi-sampling is enabled because "
                "multi-sampling runs per job after structure prediction."
            )

        # MSA (Multiple Sequence Alignment) Settings
        st.subheader("Multiple Sequence Alignment")
        enable_msa_cache = st.toggle(
            "Enable MSA Cache",
            value=True,
            help=(
                "Reuse previously generated WT/MSA data across mutant and multi-drug jobs. "
                "Recommended for screening workflows because it avoids redundant MSA server calls "
                "and significantly reduces runtime."
            ),
        )
        st.caption(
            "When enabled, BoltzOmics stores MSA from completed jobs and reuses it for related "
            "mutants/pairs. Keep this ON for mutation-heavy screens."
        )
        max_msa_seqs = st.number_input("Max MSA Sequences", min_value=1, value=8192, help="Sets the maximum number of Multiple Sequence Alignment (MSA) sequences used. Higher values improve prediction quality but increase memory usage.")
        subsample_msa = st.toggle("Subsample MSA", value=False, help="If enabled, subsamples the MSA to reduce memory usage and potentially increase prediction diversity, at the cost of slightly reduced accuracy.")
        num_subsampled_msa = 1024
        if subsample_msa:
            num_subsampled_msa = st.number_input("Number of Subsampled MSA Sequences", min_value=1, value=1024, help="Sets the number of MSA sequences to subsample. Lower values increase diversity but may reduce accuracy (recommended: 512-2048).")

        # Error handling and retry settings
        st.subheader(":material/error_outline: Error Handling")
        enable_retries = st.toggle("Enable Retries", value=True, help="Automatically retry failed predictions with exponential backoff")
        max_retry_attempts = st.slider("Max Retry Attempts", min_value=1, max_value=5, value=2, help="Maximum number of retry attempts per prediction", disabled=not enable_retries)
        retry_delay_base = st.slider("Base Retry Delay (seconds)", min_value=1, max_value=30, value=5, help="Base delay between retry attempts (doubles with each retry)", disabled=not enable_retries)
        prediction_timeout_minutes = st.slider("Job Timeout (minutes)", min_value=5, max_value=120, value=20, help="Maximum runtime per Boltz prediction before forcing a timeout.")
    
    # Page header
    with st.expander("Help", expanded=True, icon=":material/help:"):
        st.markdown("""
        **Drug Screening Workflow**

        Screen multiple drug candidates against protein targets using AI-powered structure prediction.

        **Quick Start:**
        1. **Project**: Select existing or create new project to organize results
        2. **Inputs**: Enter protein sequences (FASTA/UniProt ID/gene name) and drug SMILES
        3. **Configure**: Adjust Boltz-2 prediction parameters (optional - defaults work well)
        4. **Run**: Click "Run Drug Screening" and monitor progress

        **Input Formats:**
        - Proteins: FASTA format, UniProt accessions (e.g., P20813), or gene names (e.g., CYP2B6)
        - Drugs: SMILES strings in FASTA format (one per line with >identifier)
        - Mutations: Automatically discovered from databases or specify manually (e.g., I328T)

        **Results:** Binding predictions, confidence scores, 3D structures, and mutation impact analysis
        """)
    
    # Project management section
    with st.container(border=True):
        st.subheader(":material/folder: Project")
        
        # Get existing projects
        existing_projects = get_project_list(RESULTS_DIR)
        # Always add the option to create a new project as the first option
        project_options = ["＋ Add New Project"] + existing_projects.copy()
        # Use multiselect for project selection
        selected_projects = st.multiselect(
            "Select project(s)",
            options=project_options,
            max_selections=1,
            help="Select an existing project or choose '＋ Add New Project' to create a new one",
            key="project_selector",
        )
        
        # Handle project selection
        project_name = None
        new_project_name = None
        
        if selected_projects:
            selected = selected_projects[0]
            if selected == "＋ Add New Project":
                # Allow user to input new project name
                new_project_name = st.text_input(
                    "New project name",
                    placeholder="Enter a new project name (e.g., my_drug_screen)",
                    help="This will be the folder name where all results are stored"
                )
                if new_project_name:
                    project_name = new_project_name
            else:
                # Existing project selected
                project_name = selected
                
                # Load project data automatically
                project_data = load_project_data(project_name, RESULTS_DIR)
                if project_data:
                    # Handle case where project_data might be a list (old format)
                    if isinstance(project_data, list):
                        # Convert old format to new format
                        project_data = {
                            'project_name': project_name,
                            'created_date': datetime.now().isoformat(),
                            'last_updated': datetime.now().isoformat(),
                            'total_results': len(project_data),
                            'successful_results': len([r for r in project_data if isinstance(r, dict) and r.get("status") == "Success"]),
                            'failed_results': len([r for r in project_data if isinstance(r, dict) and r.get("status") != "Success"]),
                            'computation_time_seconds': None,  # Not available in old format
                            'results': project_data
                        }
                    
                    # Store in session state
                    st.session_state.loaded_project_data = project_data
                    st.session_state.loaded_project_name = project_name
                    # Auto-detect structure_only mode based on results
                    results_list = project_data['results'] if isinstance(project_data, dict) and 'results' in project_data else project_data
                    if results_list and all((not r.get('drug_name')) for r in results_list if isinstance(r, dict)):
                        st.session_state['structure_only'] = True
                    else:
                        st.session_state['structure_only'] = False
                    # Show simple project info
                    st.info("While the project is loaded, running additional predictions will append them to this project.")
                    
                    # Store results in session state but don't display them yet
                    if isinstance(project_data, dict) and 'results' in project_data and project_data['results']:
                        # Deduplicate results when loading from project data
                        st.session_state.screening_results = deduplicate_results(project_data['results'])
                    elif isinstance(project_data, list) and project_data:
                        # Deduplicate results when loading from project data (old format)
                        st.session_state.screening_results = deduplicate_results(project_data)
                else:
                    st.error(f"Failed to load project: {project_name}")
        if not project_name:
            project_name = st.session_state.get("loaded_project_name")
        elif project_name != st.session_state.get("loaded_project_name"):
            st.session_state.loaded_project_name = project_name
        
        # Project management buttons
        col1, col2, col3, col4, col5 = st.columns([2, 2, 2, 2, 2])
        with col1:
            deletion_mode = st.toggle("Enable Deletion Mode", value=False, help="Toggle to enable project deletion functionality")
        with col2:
            if st.button("Delete selected project", icon=":material/delete:", type="tertiary", help="Delete selected project", disabled=not deletion_mode):
                if project_name and project_name in existing_projects:
                    if delete_project(project_name, RESULTS_DIR):
                        st.success(f"Deleted project: {project_name}")
                        # Clear session state if this was the loaded project
                        if hasattr(st.session_state, 'loaded_project_name') and st.session_state.loaded_project_name == project_name:
                            del st.session_state.loaded_project_name
                            del st.session_state.loaded_project_data
                        st.rerun()
                    else:
                        st.error(f"Failed to delete project: {project_name}")
                else:
                    st.error("Please select an existing project to delete")
        
        with col3:
            if st.button("Clear loaded project", icon=":material/clear:", type="tertiary", help="Clear loaded project from memory"):
                if hasattr(st.session_state, 'loaded_project_name'):
                    del st.session_state.loaded_project_name
                if hasattr(st.session_state, 'loaded_project_data'):
                    del st.session_state.loaded_project_data
                if hasattr(st.session_state, 'screening_results'):
                    del st.session_state.screening_results
                st.success("Cleared loaded project")
                st.rerun()
        
        with col4:
            if st.button("Refresh project list", icon=":material/refresh:", type="tertiary", help="Refresh the list of existing projects"):
                st.rerun()
        
        with col5:
            if st.button("Force update project", icon=":material/sync:", type="tertiary", help="Force update project metadata from boltz results folders"):
                if project_name and project_name in existing_projects:
                    force_update_project_from_boltz_results(project_name, RESULTS_DIR)
                    # Reload project data to reflect changes
                    project_data = load_project_data(project_name, RESULTS_DIR)
                    if project_data and isinstance(project_data, dict) and 'results' in project_data:
                        st.session_state.screening_results = deduplicate_results(project_data['results'])
                    st.rerun()
                else:
                    st.error("Please select an existing project to force update")
    
    # Main content area with containers (shown for both new projects and adding to existing projects)
    with st.container(border=True):
        st.subheader(":material/input: Input Data")
        # Add segmented control for input mode and structural-only toggle
        col_input_mode, col_struct_only = st.columns(2, vertical_alignment="bottom")
        with col_input_mode:
            input_mode = st.segmented_control(
                "Input mode",
                options=["Multi-Protein Mode", "Mutation Mode"],
                default=["Multi-Protein Mode"],
                help="Multi-Protein Mode: Screen different proteins against drugs. Mutation Mode: Analyze mutations of a single wild-type protein against drugs."
            )
        with col_struct_only:
            structure_only = st.toggle(
                "Structural prediction only (no ligand)",
                key="structure_only",
                help="If enabled, skip ligand input and affinity prediction. Only protein structure will be predicted."
            )
        mutation_chain_starts: Dict[str, int] = {}
        
        if input_mode == "Multi-Protein Mode":
            # Multi-Protein mode - existing FASTA format
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader(":material/immunology: Protein Sequences (FASTA Format)")
                protein_input = st.text_area(
                    "Protein Sequences",
                    height=300,
                    placeholder=">Protein_A\nMVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ\n\n>Protein_B\nMKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG",
                    key="protein_input_normal",
                    label_visibility="collapsed"
                )
            
            with col2:
                st.subheader(":material/mixture_med: Drug SMILES (FASTA Format)")
                drug_input = st.text_area(
                    "Drug SMILES",
                    height=300,
                    placeholder=">Aspirin\nCC(=O)OC1=CC=CC=C1C(=O)O\n\n>Ibuprofen\nCC(C)CC1=CC=C(C=C1)C(C)C(=O)O\n\n>Paracetamol\nCC(=O)NC1=CC=C(O)C=C1",
                    key="drug_input_normal",
                    label_visibility="collapsed",
                    disabled=structure_only
                )
            
            # Parse inputs for multi-protein mode
            protein_sequences = []
            drug_smiles = []
            
            if protein_input.strip():
                protein_sequences = parse_fasta_sequences(protein_input)
            
            if not structure_only and drug_input.strip():
                drug_smiles = parse_smiles_list(drug_input)
        
        else:  # Mutation Mode
            # Mutation mode - wild-type sequence + mutations
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader(":material/immunology: Wild-Type Protein Sequence")
                wt_protein_input = st.text_area(
                    "Wild-Type Protein Sequence",
                    height=200,
                    placeholder="MVTPEGNVSLVDESLLVGVTDEDRAVRSAHQFYERLIGLWAPAVMEAAHELGVFAALAEAPADSGELARRLDCDARAMRVLLDALYAYDVIDRIHDTNGFRYLLSAEARECLLPGTLFSLVGKFMHDINVAWPAWRNLAEVVRHGARDTSGAESPNGIAQEDYESLVGGINFWAPPIVTTLSRKLRASGRSGDATASVLDVGCGTGLYSQLLLREFPRWTATGLDVERIATLANAQALRLGVEERFATRAGDFWRGGWGTGYDLVLFANIFHLQTPASAVRLMRHAAACLAPDGLVAVVDQIVDADREPKTPQDRFALLFAASMTNTGGGDAYTFQEYEEWFTAAGLQRIETLDTPMHRILLARRATEPSAVPEGQASENLYFQ",
                    key="wt_protein_input",
                    label_visibility="collapsed",
                    help="Enter the wild-type protein sequence (single chain or multi-chain with : separator)"
                )
                
                # Residue numbering section
                # Parse chains to get chain IDs
                chain_starts = {}
                if wt_protein_input.strip():
                    is_valid, error_msg, chains_dict, upper_seq = validate_protein_sequence(wt_protein_input.strip())
                    if is_valid and chains_dict:
                        # Create input fields for each chain
                        for chain_id in sorted(chains_dict.keys()):
                            start_number = st.number_input(
                                f"Chain {chain_id} starting residue number",
                                min_value=1,
                                value=1,
                                key=f"chain_start_{chain_id}",
                                help=f"Enter the residue number for the first residue in chain {chain_id}"
                            )
                            chain_starts[chain_id] = start_number
                        mutation_chain_starts = dict(chain_starts)
                    else:
                        # Single chain
                        start_number = st.number_input(
                            "Starting residue number",
                            min_value=1,
                            value=1,
                            key="chain_start_single",
                            help="Enter the residue number for the first residue"
                        )
                        chain_starts['A'] = start_number
                        mutation_chain_starts = dict(chain_starts)
                
                # Mutations section
                mutations_input = st.text_area(
                    "Mutations",
                    height=100,
                    placeholder="A102S,G99E,K23R/H50L",
                    key="mutations_input",
                    help="Enter mutations in format: <Wild Type Residue><Residue Number><New Residue>. Separate multiple mutants with commas. Combine multiple mutations in one mutant with slashes."
                )
                
                # Show mutation format help
                with st.popover("Mutation Format Help"):
                    st.markdown("""
                    **Mutation Format:**
                    - Single mutation: `A102S` (Wild-type residue Ala at position 102 mutated to Ser)
                    - Multiple mutations in one mutant: `K23R/H50L` (K23R AND H50L in one mutant)
                    - Multiple mutants: `A102S,G99E,K23R/H50L` (3 separate mutants)
                    
                    **Multi-Chain Proteins:**
                    - Residue numbers are automatically mapped to the correct chain based on chain start positions
                    - Example: If Chain A starts at residue 1 and Chain B starts at residue 699, then:
                      - Residue 102 → Chain A (position 101 in chain)
                      - Residue 750 → Chain B (position 51 in chain)
                    - The system will automatically detect which chain each mutation belongs to
                    """)
                
                # Add mutation verification display
                if mutations_input.strip() and wt_protein_input.strip():
                    # Initialize chain_starts for verification
                    chain_starts_verify = {}
                    is_valid_wt, _, chains_dict_wt, upper_seq_wt = validate_protein_sequence(wt_protein_input.strip())
                    if is_valid_wt and chains_dict_wt:
                        # Multi-chain protein
                        for chain_id in sorted(chains_dict_wt.keys()):
                            chain_start_key = f"chain_start_{chain_id}"
                            if chain_start_key in st.session_state:
                                chain_starts_verify[chain_id] = st.session_state[chain_start_key]
                            else:
                                chain_starts_verify[chain_id] = 1  # Default
                    else:
                        # Single chain protein
                        if "chain_start_single" in st.session_state:
                            chain_starts_verify['A'] = st.session_state["chain_start_single"]
                        else:
                            chain_starts_verify['A'] = 1  # Default
                    
                    # Create a unique key for verification that changes when chain starts change
                    chain_starts_key = str(sorted(chain_starts_verify.items()))
                    
                    # Verify mutations
                    verification_results = verify_mutations_with_wt_residues(upper_seq_wt, mutations_input, chain_starts_verify, chains_dict_wt)
                    
                    if verification_results:
                        st.markdown("##### :material/verified: Mutation Verification")
                        
                        # Create verification data for display
                        verification_data = []
                        # Determine which mutant each verification result belongs to
                        mutant_strings = [s.strip() for s in st.session_state.get("mutations_input", "").split(',')]
                        mutant_idx_map = {}
                        for idx, mutant_str in enumerate(mutant_strings):
                            # Each mutant_str can have multiple mutations (split by /)
                            mutation_parts = [s.strip() for s in mutant_str.split('/')]
                            for part in mutation_parts:
                                if not part or len(part) < 5:
                                    continue
                                wt_residue = part[0].upper()
                                residue_part = part[1:-1]
                                new_residue = part[-1].upper()
                                try:
                                    residue_number = int(residue_part)
                                    if residue_number > 0:
                                        # Use a tuple to identify the mutation
                                        mutant_idx_map[(wt_residue, residue_number, new_residue)] = idx + 1  # 1-based
                                except ValueError:
                                    continue
                        for result in verification_results:
                            is_valid, chain_id, residue_number, wt_residue, actual_residue, new_residue, context = result
                            # Find the mutant index for this mutation
                            mutant_number = mutant_idx_map.get((wt_residue, residue_number, new_residue), 1)
                            verification_data.append({
                                "#": mutant_number,
                                "Mutation": f"{wt_residue}{residue_number}{new_residue}",
                                "Chain": chain_id,
                                "Status": "Valid" if is_valid else "Invalid",
                                "Expected": wt_residue,
                                "Found": context,  # Use context instead of just the residue
                                "New": new_residue
                            })
                        
                        # Display verification table
                        verification_df = pd.DataFrame(verification_data)
                        
                        def highlight_verification_status(val):
                            if "Valid" in val:
                                return "background-color: #d4edda; color: #155724"
                            elif "Invalid" in val:
                                return "background-color: #f8d7da; color: #721c24"
                            return ""
                        
                        styled_verification_df = verification_df.style.map(highlight_verification_status, subset=['Status'])
                        
                        st.dataframe(
                            styled_verification_df,
                            use_container_width=True,
                            hide_index=True,
                            column_config={
                                "Mutation": st.column_config.TextColumn("Mutation", width="medium"),
                                "Chain": st.column_config.TextColumn("Chain", width="small"),
                                "Status": st.column_config.TextColumn("Status", width="small"),
                                "Expected": st.column_config.TextColumn("Expected WT", width="small"),
                                "Found": st.column_config.TextColumn("Found in Sequence", width="medium", help="Shows nearby residues for context"),
                                "New": st.column_config.TextColumn("New Residue", width="small")
                            },
                            key=f"verification_table_{chain_starts_key}"  # Unique key that changes with chain starts
                        )
                        
                        # Show summary
                        valid_count = len([r for r in verification_results if r[0]])
                        invalid_count = len([r for r in verification_results if not r[0]])
                        
                        if invalid_count > 0:
                            st.warning(f"{invalid_count} mutation(s) have mismatched wild-type residues. Please check your input.")
                        else:
                            st.success(f"All {valid_count} mutation(s) are valid!")
                    
                    # Store verification results in session state for use in parsed input
                    st.session_state.mutation_verification_results = verification_results
            
            with col2:
                st.subheader(":material/mixture_med: Drug SMILES (FASTA Format)")
                drug_input = st.text_area(
                    "Drug SMILES",
                    height=300,
                    placeholder=">Aspirin\nCC(=O)OC1=CC=CC=C1C(=O)O\n\n>Ibuprofen\nCC(C)CC1=CC=C(C=C1)C(C)C(=O)O\n\n>Paracetamol\nCC(=O)NC1=CC=C(O)C=C1",
                    key="drug_input_mutation",
                    label_visibility="collapsed",
                    disabled=structure_only
                )
            
            # Parse inputs for mutation mode
            protein_sequences = []
            drug_smiles = []
            
            if wt_protein_input.strip() and mutations_input.strip():
                # Parse wild-type sequence
                is_valid, error_msg, chains_dict, upper_seq = validate_protein_sequence(wt_protein_input.strip())
                
                if is_valid:
                    # Initialize chain_starts for mutation application
                    chain_starts = {}
                    if chains_dict:
                        # Multi-chain protein
                        for chain_id in sorted(chains_dict.keys()):
                            # Try to get the chain start from session state
                            chain_start_key = f"chain_start_{chain_id}"
                            if chain_start_key in st.session_state:
                                chain_starts[chain_id] = st.session_state[chain_start_key]
                            else:
                                chain_starts[chain_id] = 1  # Default
                    else:
                        # Single chain protein
                        if "chain_start_single" in st.session_state:
                            chain_starts['A'] = st.session_state["chain_start_single"]
                        else:
                            chain_starts['A'] = 1  # Default
                    mutation_chain_starts = dict(chain_starts)
                    
                    # Parse mutations
                    mutation_lists = parse_mutations(mutations_input)
                    
                    # Create wild-type entry
                    protein_sequences.append(("WT", upper_seq))
                    
                    # Create mutant entries
                    # Split by comma to get individual mutants
                    mutant_strings = [s.strip() for s in mutations_input.split(',')]
                    
                    for i, mutant_str in enumerate(mutant_strings):
                        if not mutant_str:
                            continue
                            
                        # Generate mutant name from this specific mutant string
                        mutant_name = generate_mutant_name_from_text(mutant_str)
                        
                        # Get the mutations for this specific mutant
                        if i < len(mutation_lists):
                            mutations = mutation_lists[i]
                            # Apply mutations to sequence
                            mutated_seq = apply_mutations_to_sequence(upper_seq, mutations, chain_starts, chains_dict)
                            protein_sequences.append((mutant_name, mutated_seq))
                else:
                    st.error(f"Invalid wild-type sequence: {error_msg}")
            
            if not structure_only and drug_input.strip():
                drug_smiles = parse_smiles_list(drug_input)

    # Advanced options tabs container
    with st.container(border=True):
        st.subheader(":material/tune: Advanced Options")
        mutation_steering_enabled = False
        mutation_steering_window = 3
        mutation_steering_max_distance = 6.0
        mutation_steering_apply_to_wt = False
        mutation_steering_enable_potentials = True
        method_prior_label = "None"
        method_prior_value: Optional[str] = None
        compute_mutation_local_consistency = False

        tab_mutation, tab1, tab2, tab3, tab4, tab5 = st.tabs([
            ":material/biotech: Mutation Discovery",
            ":material/hexagon: Co-factors",
            ":material/description: Structural Template",
            ":material/donut_large: Binding Pocket",
            ":material/science: Post-translational Modifications",
            ":material/filter_alt: Protein-Drug Pairing"
        ])

        # Tab 0: Mutation Discovery
        with tab_mutation:
            display_mutation_discovery_section()
            with st.expander("Mutation-Conditioned Pocket Steering", expanded=False):
                st.caption(
                    "Adds mutation-neighborhood pocket contacts automatically using your residue numbering "
                    "and chain start values, then optionally enables Boltz inference potentials."
                )
                mutation_steering_enabled = st.toggle(
                    "Enable Mutation Steering",
                    value=False,
                    disabled=structure_only or input_mode != "Mutation Mode",
                    help=(
                        "For each mutant, automatically build pocket contacts around mutated residues. "
                        "This helps keep ligand placement consistent near mutation sites."
                    ),
                )
                mutation_steering_window = st.slider(
                    "Neighborhood Window (± residues)",
                    min_value=0,
                    max_value=12,
                    value=3,
                    disabled=not mutation_steering_enabled,
                    help="How many residues around each mutation to include as pocket contacts.",
                )
                mutation_steering_max_distance = st.number_input(
                    "Steering Max Distance (Å)",
                    min_value=4.0,
                    max_value=20.0,
                    value=6.0,
                    step=0.5,
                    disabled=not mutation_steering_enabled,
                    help="Distance used for automatically generated steering contacts.",
                )
                mutation_steering_apply_to_wt = st.toggle(
                    "Apply Steering to WT",
                    value=False,
                    disabled=not mutation_steering_enabled,
                    help="If OFF, steering is applied to mutant proteins only.",
                )
                mutation_steering_enable_potentials = st.toggle(
                    "Enable Boltz Potentials (`--use_potentials`)",
                    value=True,
                    disabled=not mutation_steering_enabled,
                    help=(
                        "Boltz inference potentials add soft geometric guidance from constraints "
                        "(such as pocket contacts) during sampling."
                    ),
                )
                method_prior_label = st.selectbox(
                    "Method Prior (`--method`, optional)",
                    options=["None", "electron microscopy", "x-ray diffraction", "other"],
                    index=0,
                    disabled=not mutation_steering_enabled,
                    help=(
                        "Optional Boltz structural-method prior. Leave as None unless you are testing "
                        "target-specific ablations."
                    ),
                )
                method_prior_value = None if method_prior_label == "None" else method_prior_label
                if structure_only:
                    st.info("Mutation steering is disabled in structure-only mode.")
                elif input_mode != "Mutation Mode":
                    st.info("Switch to Mutation Mode to enable mutation-conditioned steering.")
                compute_mutation_local_consistency = st.toggle(
                    "Compute Mutation-Local Consistency Score",
                    value=False,
                    disabled=structure_only or input_mode != "Mutation Mode",
                    help=(
                        "Compares each mutant pose against WT (same drug) using interaction fingerprints "
                        "and disruption metrics, then reports a 0-1 local consistency score."
                    ),
                )

        # Tab 1: Co-factors
        with tab1:
            # Show current cofactor info if available
            current_cofactors = st.session_state.get('cofactor_info', [])
            if current_cofactors and len(current_cofactors) > 0:
                st.info(f"Current cofactors: {len(current_cofactors)} configured")
                for i, cofactor in enumerate(current_cofactors):
                    if cofactor.get('smiles'):
                        st.write(f"  Co-factor {i+1}: SMILES - {cofactor['smiles']}")
                    elif cofactor.get('ccd'):
                        st.write(f"  Co-factor {i+1}: CCD Code - {cofactor['ccd']}")

            # Number of cofactors dropdown
            num_cofactors = st.selectbox(
                "Number of co-factors to add",
                options=[0, 1, 2, 3, 4],
                index=len(current_cofactors) if current_cofactors else 0,
                help="Select how many co-factors you want to include in all predictions"
            )

            cofactors_list = []

            if num_cofactors > 0:
                for i in range(num_cofactors):
                    chain_id = chr(ord('T') + i)  # T, U, V, W
                    # Get current cofactor data if available
                    current_cofactor = current_cofactors[i] if i < len(current_cofactors) else None

                    # Add segmented control for co-factor input method
                    col1, col2 = st.columns([1, 5])
                    with col1:
                        cofactor_input_method = st.segmented_control(
                            f"Co-factor {i+1} input method",
                            options=["SMILES", "CCD Code"],
                            default=["CCD Code"],
                            help=f"Choose how to input co-factor {i+1}",
                            key=f"screening_cofactor_method_{i}"
                        )

                    with col2:
                        if cofactor_input_method == "SMILES":
                            current_smiles = current_cofactor.get('smiles', '') if current_cofactor else ''
                            cofactor_smiles = st.text_input(
                                f"Co-factor {i+1} SMILES",
                                value=current_smiles,
                                placeholder="e.g., CC(=N)Oc1ccccc1C(=O)S",
                                help=f"Enter SMILES string for co-factor {i+1}",
                                key=f"screening_cofactor_smiles_{i}"
                            )

                            # Real-time validation for SMILES
                            if cofactor_smiles.strip():
                                if validate_smiles(cofactor_smiles.strip()):
                                    st.success(f":material/task_alt: Valid SMILES: {cofactor_smiles.strip()}")
                                else:
                                    st.error(f":material/error: Invalid SMILES format: {cofactor_smiles.strip()}")

                            cofactor_ccd = ""

                        elif cofactor_input_method == "CCD Code":
                            current_ccd = current_cofactor.get('ccd', '') if current_cofactor else ''
                            cofactor_ccd = st.text_input(
                                f"Co-factor {i+1} CCD Code",
                                value=current_ccd,
                                placeholder="e.g., HEM, ATP, NAD",
                                help=f"Enter CCD code for co-factor {i+1}",
                                key=f"screening_cofactor_ccd_{i}"
                            )

                            # Real-time validation for CCD code
                            if cofactor_ccd.strip():
                                if validate_ccd_code(cofactor_ccd.strip()):
                                    st.success(f":material/task_alt: Valid CCD code: {cofactor_ccd.strip().upper()}")
                                else:
                                    st.error(f":material/error: Invalid CCD code format: {cofactor_ccd.strip()}")

                            cofactor_smiles = ""

                        # Validate and prepare co-factor info for this cofactor
                        cofactor_info = None
                        if cofactor_input_method == "SMILES" and cofactor_smiles.strip():
                            if validate_smiles(cofactor_smiles.strip()):
                                cofactor_info = {'smiles': cofactor_smiles.strip()}
                        elif cofactor_input_method == "CCD Code" and cofactor_ccd.strip():
                            if validate_ccd_code(cofactor_ccd.strip()):
                                cofactor_info = {'ccd': cofactor_ccd.strip().upper()}

                        if cofactor_info:
                            cofactors_list.append(cofactor_info)

                # Show help popovers
                col1, col2 = st.columns(2)
                with col1:
                    with st.popover("Common CCD Codes"):
                        st.markdown("""
                        **Common Co-factor CCD Codes:**
                        - **HEM**: Heme (protoporphyrin IX)
                        - **ATP**: Adenosine triphosphate
                        - **ADP**: Adenosine diphosphate
                        - **AMP**: Adenosine monophosphate
                        - **NAD**: Nicotinamide adenine dinucleotide
                        - **NADP**: Nicotinamide adenine dinucleotide phosphate
                        - **FAD**: Flavin adenine dinucleotide
                        - **FMN**: Flavin mononucleotide
                        - **COA**: Coenzyme A
                        - **PLP**: Pyridoxal 5'-phosphate
                        - **THF**: Tetrahydrofolate
                        - **B12**: Vitamin B12 (cobalamin)
                        """)

            # Store cofactor info in session state
            st.session_state.cofactor_info = cofactors_list

            # Set cofactor_info to the list of cofactors for use in screening prediction
            cofactor_info = cofactors_list if cofactors_list else None

        # Tab 2: Structural Template
        with tab2:
            template_cif_path = None
            col, _ = st.columns([1, 1.5])
            with col:
                uploaded_cif = st.file_uploader("Upload structural template file", type=["cif"], help="Upload a .cif protein structural file to use as a template for structure prediction.")
                if uploaded_cif is not None:
                    if project_name:
                        project_dir = os.path.join(RESULTS_DIR, project_name)
                        os.makedirs(project_dir, exist_ok=True)
                        cif_filename = f"template_{datetime.now().strftime('%Y%m%d_%H%M%S')}.cif"
                        cif_path = os.path.join(project_dir, cif_filename)
                        with open(cif_path, "wb") as f:
                            f.write(uploaded_cif.read())
                        template_cif_path = os.path.abspath(cif_path)
                        st.session_state["template_cif_path"] = template_cif_path
                    else:
                        st.warning("Please select or create a project before uploading a template file.")
                else:
                    # If nothing uploaded, clear session state for template path
                    st.session_state["template_cif_path"] = None
                # For downstream use
                template_cif_path = st.session_state.get("template_cif_path")

        # Tab 3: Binding Pocket
        with tab3:
            # Binding pocket constraints section
            binding_pocket_constraints = None

            # Get the first protein sequence for binding pocket analysis if available
            protein_sequence_for_pocket = None
            if protein_sequences:
                # For mutation mode, use the wild-type sequence
                if input_mode == "Mutation Mode":
                    # Find the WT sequence
                    for name, seq in protein_sequences:
                        if name == "WT":
                            protein_sequence_for_pocket = seq
                            break
                else:
                    # For multi-protein mode, use the first protein sequence
                    protein_sequence_for_pocket = protein_sequences[0][1]  # Get sequence from first protein

            # Call display_binding_pocket_section with the protein sequence (can be None)
            binding_pocket_constraints = display_binding_pocket_section(protein_sequence_for_pocket)

            # Store binding pocket constraints in session state
            if binding_pocket_constraints:
                st.session_state.binding_pocket_constraints = binding_pocket_constraints

        # Tab 4: Post-translational Modifications
        with tab4:
            # Call display_ptm_section with the protein sequence (can be None)
            protein_sequence_for_ptm = None
            if protein_sequences:
                # For mutation mode, use the wild-type sequence
                if input_mode == "Mutation Mode":
                    # Find the WT sequence
                    for name, seq in protein_sequences:
                        if name == "WT":
                            protein_sequence_for_ptm = seq
                            break
                else:
                    # For multi-protein mode, use the first protein sequence
                    protein_sequence_for_ptm = protein_sequences[0][1]  # Get sequence from first protein

            ptm_modifications = display_ptm_section(protein_sequence_for_ptm, None, protein_sequences)

            # Store PTM modifications in session state
            if ptm_modifications:
                st.session_state.ptm_modifications = ptm_modifications

        # Tab 5: Protein-Drug Pairing
        with tab5:
            # Call display_protein_drug_filter_section
            protein_drug_filter = display_protein_drug_filter_section(protein_sequences, drug_smiles)

            # Store protein-drug filter in session state
            if protein_drug_filter:
                st.session_state.protein_drug_filter = protein_drug_filter

    mutation_steering_config = {
        "enabled": bool(mutation_steering_enabled),
        "window": int(mutation_steering_window),
        "max_distance": float(mutation_steering_max_distance),
        "apply_to_wt": bool(mutation_steering_apply_to_wt),
        "mutation_mode_active": input_mode == "Mutation Mode",
        "chain_starts": dict(mutation_chain_starts),
    }
    st.session_state["mutation_steering_config"] = mutation_steering_config
    st.session_state["compute_mutation_local_consistency"] = bool(compute_mutation_local_consistency)
    use_potentials = bool(mutation_steering_enabled and mutation_steering_enable_potentials)
    method = method_prior_value if mutation_steering_enabled else None
    
    # Display parsed data
    if protein_sequences or drug_smiles:
        with st.container(border=True):
            st.subheader(":material/preview: Parsed Input")
            
            col1, col2 = st.columns(2)
            
            with col1:
                if protein_sequences:
                    st.write(f"**:material/immunology: Proteins ({len(protein_sequences)}):**")
                    
                    # Create protein data for table
                    protein_data = []
                    for name, seq in protein_sequences:
                        # Validate protein sequence using chain parsing logic
                        is_valid, error_msg, chains_dict, upper_seq = validate_protein_sequence(seq)
                        
                        # For mutation mode, check verification results
                        if input_mode == "Mutation Mode" and name != "WT":
                            # Get verification results from session state
                            verification_results = st.session_state.get("mutation_verification_results", [])
                            
                            # Check if this mutant has any invalid mutations
                            has_invalid_mutations = False
                            if verification_results:
                                # Find mutations for this specific mutant
                                mutant_strings = [s.strip() for s in st.session_state.get("mutations_input", "").split(',')]
                                for i, mutant_str in enumerate(mutant_strings):
                                    if not mutant_str:
                                        continue
                                    
                                    mutant_name_check = generate_mutant_name_from_text(mutant_str)
                                    if mutant_name_check == name:
                                        # Parse the specific mutations for this mutant
                                        mutation_parts = [s.strip() for s in mutant_str.split('/')]
                                        mutant_mutations = []
                                        
                                        for part in mutation_parts:
                                            if not part or len(part) < 5:
                                                continue
                                            # Extract wild-type residue, residue number, and new residue
                                            wt_residue = part[0].upper()
                                            residue_part = part[1:-1]
                                            new_residue = part[-1].upper()
                                            
                                            try:
                                                residue_number = int(residue_part)
                                                if residue_number > 0:
                                                    mutant_mutations.append(f"{wt_residue}{residue_number}{new_residue}")
                                            except ValueError:
                                                continue
                                        
                                        # Check if any of this mutant's mutations are invalid
                                        for mutation_str in mutant_mutations:
                                            for result in verification_results:
                                                is_valid, chain_id, residue_number, wt_residue, actual_residue, new_residue, context = result
                                                expected_mutation = f"{wt_residue}{residue_number}{new_residue}"
                                                if expected_mutation == mutation_str and not is_valid:
                                                    has_invalid_mutations = True
                                                    break
                                            if has_invalid_mutations:
                                                break
                                        break
                            
                            if has_invalid_mutations:
                                status = "Invalid - Verification Failed"
                                length = f"{len(seq)} aa"
                                chains = "Single chain" if not chains_dict else ", ".join([f"Chain {chain_id}: {len(chain_seq)} aa" for chain_id, chain_seq in chains_dict.items()])
                            elif is_valid:
                                status = "Valid"
                                length = f"{len(seq)} aa"
                                if chains_dict:
                                    chain_info = ", ".join([f"Chain {chain_id}: {len(chain_seq)} aa" for chain_id, chain_seq in chains_dict.items()])
                                    chains = chain_info
                                else:
                                    chains = "Single chain"
                            else:
                                status = "Invalid"
                                length = "N/A"
                                chains = error_msg
                        else:
                            # Multi-protein mode or WT sequence
                            if is_valid:
                                status = "Valid"
                                length = f"{len(seq)} aa"
                                if chains_dict:
                                    chain_info = ", ".join([f"Chain {chain_id}: {len(chain_seq)} aa" for chain_id, chain_seq in chains_dict.items()])
                                    chains = chain_info
                                else:
                                    chains = "Single chain"
                            else:
                                status = "Invalid"
                                length = "N/A"
                                chains = error_msg
                        
                        protein_data.append({
                            "Name": name,
                            "Status": status,
                            "Length": length,
                            "Chains": chains,
                            "Sequence": seq[:50] + "..." if len(seq) > 50 else seq
                        })
                        
                        # For mutation mode, highlight mutations in sequence
                        if input_mode == "Mutation Mode" and name != "WT":
                            # Get mutations input from session state
                            mutations_input_display = st.session_state.get("mutations_input", "")
                            if mutations_input_display:
                                # Parse mutations for this specific mutant
                                mutation_lists = parse_mutations(mutations_input_display)
                                
                                # Find which mutant this corresponds to
                                mutant_strings = [s.strip() for s in mutations_input_display.split(',')]
                                for i, mutant_str in enumerate(mutant_strings):
                                    if not mutant_str:
                                        continue
                                        
                                    # Generate mutant name from this specific mutant string
                                    mutant_name_check = generate_mutant_name_from_text(mutant_str)
                                    
                                    if mutant_name_check == name and i < len(mutation_lists):
                                        mutations = mutation_lists[i]
                                        
                                        # Get chain starts from the current session
                                        chain_starts = {}
                                        wt_protein_input_display = st.session_state.get("wt_protein_input", "")
                                        if wt_protein_input_display.strip():
                                            is_valid_wt, _, chains_dict_wt, _ = validate_protein_sequence(wt_protein_input_display.strip())
                                            if is_valid_wt and chains_dict_wt:
                                                for chain_id in sorted(chains_dict_wt.keys()):
                                                    # Try to get the chain start from session state
                                                    chain_start_key = f"chain_start_{chain_id}"
                                                    if chain_start_key in st.session_state:
                                                        chain_starts[chain_id] = st.session_state[chain_start_key]
                                                    else:
                                                        chain_starts[chain_id] = 1  # Default
                                            else:
                                                # Single chain
                                                if "chain_start_single" in st.session_state:
                                                    chain_starts['A'] = st.session_state["chain_start_single"]
                                                else:
                                                    chain_starts['A'] = 1
                                        
                                        # Create a simple text representation with mutation indicators
                                        seq_list = list(seq)
                                        mutation_positions = set()
                                        
                                        for chain_id, residue_number, new_residue in mutations:
                                            # For multi-chain proteins, we need to handle both heteromultimers and homomultimers
                                            # In homomultimers, all chains start at the same residue number and mutations should be applied to all chains
                                            target_chains = []
                                            if chains_dict_wt and chain_id == 'A' and len(chains_dict_wt) > 1:
                                                # Find all chains that contain this residue number
                                                for cid in sorted(chains_dict_wt.keys()):
                                                    if cid in chain_starts:
                                                        chain_start = chain_starts[cid]
                                                        chain_length = len(chains_dict_wt[cid])
                                                        chain_end = chain_start + chain_length - 1
                                                        
                                                        if chain_start <= residue_number <= chain_end:
                                                            target_chains.append(cid)
                                            else:
                                                # Specific chain was specified or single chain protein
                                                target_chains = [chain_id]
                                            
                                            # Apply mutation display to all target chains
                                            for target_chain_id in target_chains:
                                                if target_chain_id in chain_starts and target_chain_id in chains_dict_wt:
                                                    # Calculate position within the specific chain
                                                    chain_start = chain_starts[target_chain_id]
                                                    chain_sequence = chains_dict_wt[target_chain_id]
                                                    position_in_chain = residue_number - chain_start
                                                    
                                                    # Check if position is valid within the chain
                                                    if 0 <= position_in_chain < len(chain_sequence):
                                                        # Calculate position in concatenated sequence
                                                        # We need to find the offset of this chain in the concatenated sequence
                                                        offset = 0
                                                        for prev_chain_id, prev_chain_seq in chains_dict_wt.items():
                                                            if prev_chain_id == target_chain_id:
                                                                break
                                                            offset += len(prev_chain_seq)
                                                        
                                                        absolute_pos = offset + position_in_chain
                                                        if 0 <= absolute_pos < len(seq_list):
                                                            mutation_positions.add(absolute_pos)
                                        
                                        # Create sequence display with mutation indicators
                                        display_parts = []
                                        for i, residue in enumerate(seq_list):
                                            if i in mutation_positions:
                                                display_parts.append(f"[{residue}]")  # Brackets for mutations
                                            else:
                                                display_parts.append(residue)
                                        
                                        # Update the sequence display in the protein data
                                        sequence_display = ''.join(display_parts)
                                        if len(sequence_display) > 50:
                                            sequence_display = sequence_display[:50] + "..."
                                        
                                        # Update the last added protein data entry
                                        if protein_data:
                                            protein_data[-1]["Sequence"] = sequence_display
                                        break
                    
                    # Display protein table with color highlighting
                    protein_df = pd.DataFrame(protein_data)
                    
                    # Apply color highlighting
                    def highlight_status(val):
                        if val == "Valid":
                            return "background-color: #d4edda; color: #155724"  # Light green background, dark green text
                        elif val == "Invalid" or val == "Invalid - Verification Failed":
                            return "background-color: #f8d7da; color: #721c24"  # Light red background, dark red text
                        else:
                            return "background-color: lightyellow"
                    
                    styled_protein_df = protein_df.style.map(highlight_status, subset=['Status'])
                    
                    st.dataframe(
                        styled_protein_df,
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "Name": st.column_config.TextColumn("Protein Name", width="medium"),
                            "Status": st.column_config.TextColumn("Status", width="small"),
                            "Length": st.column_config.TextColumn("Length", width="small"),
                            "Chains": st.column_config.TextColumn("Chain Info", width="medium"),
                            "Sequence": st.column_config.TextColumn("Sequence", width="large", help="Bracketed residues [X] indicate mutations")
                        }
                    )
                else:
                    st.write("**:material/immunology: Proteins:** None parsed")
            
            with col2:
                if drug_smiles:
                    st.write(f"**:material/mixture_med: Drugs ({len(drug_smiles)}):**")
                    
                    # Create drug data for table
                    drug_data = []
                    for name, smiles in drug_smiles:
                        # Validate SMILES
                        is_valid = validate_smiles(smiles)
                        if is_valid:
                            status = "Valid"
                        else:
                            status = "Invalid"
                        
                        # Truncate SMILES for display
                        display_smiles = smiles[:50] + "..." if len(smiles) > 50 else smiles
                        
                        drug_data.append({
                            "Name": name,
                            "Status": status,
                            "SMILES": display_smiles
                        })
                    
                    # Display drug table with color highlighting
                    drug_df = pd.DataFrame(drug_data)
                    
                    # Apply color highlighting
                    def highlight_status(val):
                        if val == "Valid":
                            return "background-color: #d4edda; color: #155724"  # Light green background, dark green text
                        elif val == "Invalid":
                            return "background-color: #f8d7da; color: #721c24"  # Light red background, dark red text
                        return ""
                    
                    styled_drug_df = drug_df.style.map(highlight_status, subset=['Status'])
                    
                    st.dataframe(
                        styled_drug_df,
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "Name": st.column_config.TextColumn("Drug Name", width="medium"),
                            "Status": st.column_config.TextColumn("Status", width="small"),
                            "SMILES": st.column_config.TextColumn("SMILES", width="large", help="Click to see full SMILES")
                        }
                    )
                else:
                    st.write("**:material/mixture_med: Drugs:** None parsed")
            
            # Display protein-drug filtering verification if enabled
            protein_drug_filter = st.session_state.get('protein_drug_filter')
            if protein_drug_filter and protein_drug_filter.get('enabled'):
                filter_pairs = protein_drug_filter.get('pairs', [])
                if filter_pairs:
                    st.markdown("---")
                    st.write("**:material/filter_alt: Protein-Drug Combinations to Evaluate:**")
                    
                    # Create verification data for filtered combinations
                    verification_data = []
                    available_proteins = {name for name, _ in protein_sequences} if protein_sequences else set()
                    available_drugs = {name for name, _ in drug_smiles} if drug_smiles else set()
                    
                    for protein_name, drug_name in filter_pairs:
                        # Check if protein and drug exist in the parsed sequences
                        protein_exists = protein_name in available_proteins
                        drug_exists = drug_name in available_drugs if not structure_only else True
                        
                        # Determine status
                        if structure_only:
                            status = "✅ Valid" if protein_exists else "❌ Protein not found"
                        else:
                            if protein_exists and drug_exists:
                                status = "✅ Valid"
                            elif not protein_exists and not drug_exists:
                                status = "❌ Both not found"
                            elif not protein_exists:
                                status = "❌ Protein not found"
                            else:
                                status = "❌ Drug not found"
                        
                        verification_data.append({
                            "Protein": protein_name,
                            "Drug": drug_name if not structure_only else "N/A (Structure only)",
                            "Status": status
                        })
                    
                    # Display verification table
                    if verification_data:
                        verification_df = pd.DataFrame(verification_data)
                        
                        # Define color styling function
                        def highlight_status(val):
                            if "✅ Valid" in val:
                                return "background-color: #d4edda; color: #155724"
                            elif "❌" in val:
                                return "background-color: #f8d7da; color: #721c24"
                            return ""
                        
                        # Apply styling
                        styled_df = verification_df.style.map(highlight_status, subset=['Status'])
                        
                        st.dataframe(
                            styled_df,
                            use_container_width=True,
                            hide_index=True,
                            column_config={
                                "Protein": st.column_config.TextColumn("Protein", width="medium"),
                                "Drug": st.column_config.TextColumn("Drug", width="medium"),
                                "Status": st.column_config.TextColumn("Status", width="medium")
                            }
                        )
                        
                        # Show summary
                        valid_count = sum(1 for item in verification_data if "✅ Valid" in item["Status"])
                        total_count = len(verification_data)
                        if valid_count == total_count:
                            st.success(f"All {valid_count} protein-drug combinations are valid and will be evaluated.")
                        elif valid_count > 0:
                            st.warning(f"{valid_count} out of {total_count} protein-drug combinations are valid and will be evaluated.")
                        else:
                            st.error("No valid protein-drug combinations found. Check your filter input.")
                elif protein_drug_filter.get('enabled'):
                    st.info("**:material/filter_alt: Protein-Drug Filter:** Enabled but no pairs specified - no combinations will be evaluated.")
            
            # Display co-factor info if provided
            if cofactor_info:
                st.write("**:material/science: Co-factors:**")
                if isinstance(cofactor_info, list):
                    for i, cofactor in enumerate(cofactor_info):
                        chain_id = chr(ord('T') + i)
                        if cofactor.get('smiles'):
                            st.write(f"- Co-factor {i+1} (Chain {chain_id}): SMILES - {cofactor['smiles']}")
                        elif cofactor.get('ccd'):
                            st.write(f"- Co-factor {i+1} (Chain {chain_id}): CCD Code - {cofactor['ccd']}")
                else:
                    # Backward compatibility for single cofactor dict
                    if cofactor_info.get('smiles'):
                        st.write(f"- SMILES: {cofactor_info['smiles']}")
                    elif cofactor_info.get('ccd'):
                        st.write(f"- CCD Code: {cofactor_info['ccd']}")
                    st.write("- Chain ID: T")
            
            # Display binding pocket constraints if provided
            binding_pocket_constraints = st.session_state.get('binding_pocket_constraints')
            if binding_pocket_constraints and binding_pocket_constraints.get('contacts'):
                contacts = binding_pocket_constraints.get('contacts', [])
                max_distance = binding_pocket_constraints.get('max_distance', 7.0)
                binder = binding_pocket_constraints.get('binder', 'X')
                # Show first contact as the pocket label
                if contacts:
                    first_contact = contacts[0]
                    pocket_label = f"{first_contact[0]}{first_contact[1]}"
                    example_residues = ', '.join([f"{c[0]}{c[1]}" for c in contacts[:5]])
                    if len(contacts) > 5:
                        example_residues += f" ... and {len(contacts) - 5} more"
                    st.write(f"**:material/donut_large: Binding pocket [<chain><resid>]:**")
                    st.write(f"- {example_residues}")
                st.write(f" - Max distance: {max_distance} Å")

            steering_cfg = st.session_state.get("mutation_steering_config", {})
            if steering_cfg.get("enabled"):
                st.write("**:material/biotech: Mutation-conditioned steering:**")
                st.write(f"- Window: ±{int(steering_cfg.get('window', 3))} residues around each mutation")
                st.write(f"- Auto-generated max distance: {float(steering_cfg.get('max_distance', 6.0)):.1f} Å")
                st.write(f"- Apply to WT: {'Yes' if steering_cfg.get('apply_to_wt') else 'No'}")
                st.write(f"- Boltz potentials (`--use_potentials`): {'Enabled' if use_potentials else 'Disabled'}")
                if method:
                    st.write(f"- Method prior (`--method`): {method}")
            if not structure_only:
                st.write("**:material/sync_alt: Affinity prediction setup:**")
                if affinity_multisampling_enabled:
                    st.write("- Multi-sampling is ON (one structure run, then multiple affinity passes).")
                    profile_text = ", ".join(
                        f"{int(p['sampling_steps_affinity'])}:{int(p['diffusion_samples_affinity'])}"
                        for p in (affinity_multisampling_profiles or [])
                    )
                    st.write(
                        f"- Settings tested (steps:diffusion): "
                        f"{profile_text or f'{int(sampling_steps_affinity)}:{int(diffusion_samples_affinity)}'}"
                    )
                    if affinity_multisampling_apply_aggregate:
                        st.write(
                            f"- Final reported affinity: {str(affinity_multisampling_aggregate_mode).replace('_', ' ')} across tested settings"
                        )
                    if affinity_multisampling_early_stop_enabled:
                        st.write("- Early stop is ON (stops when results stabilize).")
                    if affinity_multisampling_robust_outlier_filter:
                        st.write("- Outlier filtering is ON before combining settings.")
                else:
                    st.write("- Single affinity setting is ON (no sweep).")
                    st.write(
                        f"- Setting (steps:diffusion): "
                        f"{int(sampling_steps_affinity)}:{int(diffusion_samples_affinity)}"
                    )

            # Display PTM modifications if provided
            ptm_modifications = st.session_state.get('ptm_modifications')
            if ptm_modifications and ptm_modifications.get('modifications'):
                modifications = ptm_modifications.get('modifications', [])
                if modifications:
                    st.write(f"**:material/science: Post-translational modifications:**")
                    for mod in modifications:
                        protein = mod.get('protein', 'All models')
                        chain_id = mod.get('chain_id', 'A')
                        position = mod.get('position', '')
                        ccd = mod.get('ccd', '')
                        st.write(f"- {protein}, Chain {chain_id}, position {position}: {ccd}")
    
    manager = get_job_manager() if USE_SCREENING_JOB_QUEUE else None
    queue_mode_active = USE_SCREENING_JOB_QUEUE and manager is not None
    with st.container():
        filter_state = st.session_state.get('protein_drug_filter')
        total_prediction_jobs = calculate_filtered_job_count(
            protein_sequences,
            drug_smiles,
            structure_only,
            filter_state,
        ) if (protein_sequences and (drug_smiles or structure_only)) else 0

        queue_summary_for_controls = (
            manager.get_project_summary(project_name)
            if (queue_mode_active and manager and project_name)
            else None
        )

        pad_left, run_col, cancel_col, pad_right = st.columns([1, 2, 2, 1])
        with run_col:
            disabled = not project_name or not protein_sequences or (not drug_smiles and not structure_only)
            if st.button(
                "Run Drug Screening",
                icon=":material/play_circle:",
                type="primary",
                use_container_width=True,
                disabled=disabled,
            ):
                if utils is None:
                    st.error(":material/error: Utils module not available. Please ensure the main Boltzomics application is properly installed.")
                    st.stop()

                binding_pocket_constraints = st.session_state.get('binding_pocket_constraints')
                ptm_modifications = st.session_state.get('ptm_modifications')

                if total_prediction_jobs == 0:
                    if filter_state and filter_state.get('enabled'):
                        st.error("The protein-drug filter excludes all combinations. No predictions will be run.")
                    else:
                        st.error("No valid protein-drug combinations found.")
                elif queue_mode_active:
                    ensure_job_manager_executor(worker_count=queue_worker_count)
                    shared_params = {
                        "use_gpu": use_gpu,
                        "accelerator": accelerator,
                        "devices": boltz_devices,
                        "cuda_visible_devices": cuda_visible_devices,
                        "queue_gpu_devices": queue_gpu_devices if gpu_execution_mode.startswith("Multi-GPU") else None,
                        "preprocessing_threads": preprocessing_threads,
                        "override": not use_existing_results,
                        "recycling_steps": recycling_steps,
                        "sampling_steps": sampling_steps,
                        "diffusion_samples": diffusion_samples,
                        "max_parallel_samples": max_parallel_samples,
                        "step_scale": step_scale,
                        "affinity_mw_correction": affinity_mw_correction,
                        "affinity_consensus_enabled": affinity_consensus_enabled,
                        "affinity_consensus_mode": affinity_consensus_mode,
                        "affinity_consensus_weight_floor": affinity_consensus_weight_floor,
                        "affinity_consensus_entropy_alpha": affinity_consensus_entropy_alpha,
                        "external_boltz_patch_enabled": external_boltz_patch_enabled,
                        "external_boltz_patch_mode": external_boltz_patch_mode,
                        "external_boltz_patch_weight_floor": external_boltz_patch_weight_floor,
                        "external_boltz_patch_entropy_alpha": external_boltz_patch_entropy_alpha,
                        "external_boltz_patch_uncertainty_penalty": external_boltz_patch_uncertainty_penalty,
                        "external_boltz_patch_min_confidence": external_boltz_patch_min_confidence,
                        "affinity_multisampling_enabled": affinity_multisampling_enabled,
                        "affinity_multisampling_profiles": affinity_multisampling_profiles,
                        "affinity_multisampling_settings": affinity_multisampling_settings,
                        "affinity_multisampling_refinement_steps": affinity_multisampling_refinement_steps,
                        "affinity_multisampling_aggregate_mode": affinity_multisampling_aggregate_mode,
                        "affinity_multisampling_apply_aggregate": affinity_multisampling_apply_aggregate,
                        "affinity_multisampling_early_stop_enabled": affinity_multisampling_early_stop_enabled,
                        "affinity_multisampling_early_stop_min_points": affinity_multisampling_early_stop_min_points,
                        "affinity_multisampling_early_stop_delta": affinity_multisampling_early_stop_delta,
                        "affinity_multisampling_early_stop_std": affinity_multisampling_early_stop_std,
                        "affinity_multisampling_early_stop_patience": affinity_multisampling_early_stop_patience,
                        "affinity_multisampling_robust_outlier_filter": affinity_multisampling_robust_outlier_filter,
                        "affinity_multisampling_robust_outlier_zmax": affinity_multisampling_robust_outlier_zmax,
                        "affinity_multisampling_bootstrap_samples": affinity_multisampling_bootstrap_samples,
                        "use_potentials": use_potentials,
                        "method": method,
                        "max_msa_seqs": max_msa_seqs,
                        "sampling_steps_affinity": sampling_steps_affinity,
                        "diffusion_samples_affinity": diffusion_samples_affinity,
                        "enable_retries": enable_retries,
                        "max_retry_attempts": max_retry_attempts,
                        "retry_delay_base": retry_delay_base,
                        "subsample_msa": subsample_msa,
                        "num_subsampled_msa": num_subsampled_msa,
                        "template_cif_path": template_cif_path,
                        "binding_pocket_constraints": binding_pocket_constraints,
                        "mutation_steering_config": mutation_steering_config,
                        "cofactor_info": cofactor_info,
                        "ptm_modifications": ptm_modifications,
                        "prediction_timeout_seconds": prediction_timeout_minutes * 60,
                        "enable_msa_cache": enable_msa_cache,
                    }
                    jobs, cached_results, job_summary = prepare_screening_jobs(
                        protein_sequences=protein_sequences,
                        drug_smiles=drug_smiles,
                        project_name=project_name,
                        structure_only=structure_only,
                        use_existing_results=use_existing_results,
                        protein_drug_filter=filter_state,
                        shared_params=shared_params,
                        manager=manager,
                    )
                    if cached_results:
                        current_results = st.session_state.get('screening_results', [])
                        st.session_state.screening_results = deduplicate_results(current_results + cached_results)
                        st.info(f"Loaded {len(cached_results)} cached result(s).")
                    if jobs:
                        manager = get_job_manager()
                        if manager is None:
                            st.error("Job queue is unavailable. Please try again.")
                        else:
                            enqueued_jobs = manager.enqueue_jobs(jobs)
                            if enqueued_jobs:
                                st.success(f"Queued {len(enqueued_jobs)} screening job(s).")
                            else:
                                st.info("All requested combinations are already processed or queued.")
                            skipped_in_enqueue = len(jobs) - len(enqueued_jobs)
                            if skipped_in_enqueue > 0:
                                st.info(f"Skipped {skipped_in_enqueue} job(s) that were already active in the queue.")
                    elif not cached_results:
                        st.warning("No new screening jobs were scheduled.")
                    for warning_msg in job_summary.get("warnings", []):
                        st.warning(warning_msg)
                    synchronize_job_results(project_name)
                    st.session_state.loaded_project_name = project_name
                else:
                    with st.spinner("Running screening prediction..."):
                        results, computation_time = run_screening_prediction(
                            protein_sequences=protein_sequences,
                            drug_smiles=drug_smiles,
                            project_name=project_name,
                            use_gpu=use_gpu,
                            accelerator=accelerator,
                            devices=boltz_devices,
                            cuda_visible_devices=cuda_visible_devices,
                            preprocessing_threads=preprocessing_threads,
                            use_existing_results=use_existing_results,
                            recycling_steps=recycling_steps,
                            sampling_steps=sampling_steps,
                            diffusion_samples=diffusion_samples,
                            max_parallel_samples=max_parallel_samples,
                            step_scale=step_scale,
                            affinity_mw_correction=affinity_mw_correction,
                            affinity_consensus_enabled=affinity_consensus_enabled,
                            affinity_consensus_mode=affinity_consensus_mode,
                            affinity_consensus_weight_floor=affinity_consensus_weight_floor,
                            affinity_consensus_entropy_alpha=affinity_consensus_entropy_alpha,
                            external_boltz_patch_enabled=external_boltz_patch_enabled,
                            external_boltz_patch_mode=external_boltz_patch_mode,
                            external_boltz_patch_weight_floor=external_boltz_patch_weight_floor,
                            external_boltz_patch_entropy_alpha=external_boltz_patch_entropy_alpha,
                            external_boltz_patch_uncertainty_penalty=external_boltz_patch_uncertainty_penalty,
                            external_boltz_patch_min_confidence=external_boltz_patch_min_confidence,
                            affinity_multisampling_enabled=affinity_multisampling_enabled,
                            affinity_multisampling_profiles=affinity_multisampling_profiles,
                            affinity_multisampling_settings=affinity_multisampling_settings,
                            affinity_multisampling_refinement_steps=affinity_multisampling_refinement_steps,
                            affinity_multisampling_aggregate_mode=affinity_multisampling_aggregate_mode,
                            affinity_multisampling_apply_aggregate=affinity_multisampling_apply_aggregate,
                            affinity_multisampling_early_stop_enabled=affinity_multisampling_early_stop_enabled,
                            affinity_multisampling_early_stop_min_points=affinity_multisampling_early_stop_min_points,
                            affinity_multisampling_early_stop_delta=affinity_multisampling_early_stop_delta,
                            affinity_multisampling_early_stop_std=affinity_multisampling_early_stop_std,
                            affinity_multisampling_early_stop_patience=affinity_multisampling_early_stop_patience,
                            affinity_multisampling_robust_outlier_filter=affinity_multisampling_robust_outlier_filter,
                            affinity_multisampling_robust_outlier_zmax=affinity_multisampling_robust_outlier_zmax,
                            affinity_multisampling_bootstrap_samples=affinity_multisampling_bootstrap_samples,
                            use_potentials=use_potentials,
                            method=method,
                            max_msa_seqs=max_msa_seqs,
                            sampling_steps_affinity=sampling_steps_affinity,
                            diffusion_samples_affinity=diffusion_samples_affinity,
                            cofactor_info=cofactor_info,
                            binding_pocket_constraints=binding_pocket_constraints,
                            mutation_steering_config=mutation_steering_config,
                            enable_retries=enable_retries,
                            max_retry_attempts=max_retry_attempts,
                            retry_delay_base=retry_delay_base,
                            subsample_msa=subsample_msa,
                            num_subsampled_msa=num_subsampled_msa,
                            template_cif_path=template_cif_path,
                            structure_only=structure_only,
                            ptm_modifications=ptm_modifications,
                            prediction_timeout_seconds=prediction_timeout_minutes * 60,
                            enable_msa_cache=enable_msa_cache,
                            enable_batch_execution=(enable_batch_execution and not affinity_multisampling_enabled),
                        )
                    if results:
                        current_results = st.session_state.get('screening_results', [])
                        st.session_state.screening_results = deduplicate_results(current_results + results)
                    st.session_state.last_computation_time = computation_time
                    project_data = load_project_data(project_name, RESULTS_DIR)
                    if project_data:
                        if isinstance(project_data, dict) and 'results' in project_data:
                            st.session_state.screening_results = deduplicate_results(project_data['results'])
                        elif isinstance(project_data, list):
                            st.session_state.screening_results = deduplicate_results(project_data)
                        else:
                            st.session_state.screening_results = []
                    else:
                        st.session_state.screening_results = []
                    st.session_state.loaded_project_name = project_name
                    st.success(f":material/check_circle: Drug screening completed! Processed {total_prediction_jobs} combinations.")

        with cancel_col:
            has_pending = (
                queue_mode_active
                and manager
                and queue_summary_for_controls
                and queue_summary_for_controls.get("pending", 0) > 0
            )
            if st.button(
                "Cancel Pending Jobs",
                icon=":material/cancel:",
                type="tertiary",
                use_container_width=True,
                disabled=not has_pending,
                key=f"cancel_jobs_inline_{project_name}"
            ):
                cancel_info = manager.cancel_project_jobs(project_name)
                if cancel_info.get("pending_removed"):
                    st.warning(f"Cancelled {cancel_info['pending_removed']} pending job(s).")
                else:
                    st.info("No pending jobs were waiting in the queue.")
                st.rerun()

        if not project_name:
            st.info(":material/info: Please select or create a project folder.")

    if queue_mode_active and project_name:
        synchronize_job_results(project_name)
        queue_summary = render_job_queue_status(project_name)
        maybe_schedule_queue_autorefresh(queue_summary)

    # Optional post-processing: mutation-local consistency annotation.
    active_project_name = project_name or st.session_state.get("loaded_project_name")
    if (
        annotate_results_with_mutation_local_consistency is not None
        and not structure_only
        and active_project_name
        and hasattr(st.session_state, "screening_results")
        and st.session_state.screening_results
        and st.session_state.get("compute_mutation_local_consistency", False)
    ):
        cache = st.session_state.get("_mutation_local_consistency_cache", {})
        annotated_results, cache, consistency_summary = annotate_results_with_mutation_local_consistency(
            results=st.session_state.screening_results,
            project_name=active_project_name,
            enabled=True,
            cache=cache,
        )
        st.session_state._mutation_local_consistency_cache = cache
        st.session_state.screening_results = annotated_results
        if consistency_summary.get("annotated", 0) > 0:
            st.caption(
                "Mutation-local consistency annotated: "
                f"{consistency_summary['annotated']} entries "
                f"(errors: {consistency_summary.get('errors', 0)})."
            )

    # Display screening processing summary outside of columns
    if hasattr(st.session_state, 'screening_results') and st.session_state.screening_results:
        # Count existing vs new results if we have computation time
        if hasattr(st.session_state, 'last_computation_time') and st.session_state.last_computation_time:
            # This is a rough estimate - in a real scenario you'd track this more precisely
            total_results = len(st.session_state.screening_results)
            # Always use st.session_state.screening_results for summary
            st.info(f"Screening summary: {total_results} total results loaded")
        else:
            # Show basic summary when computation time is not available
            total_results = len(st.session_state.screening_results)
            successful = len([r for r in st.session_state.screening_results if r.get("status") == "Success"])
            st.info(f"Screening summary: {total_results} total results ({successful} successful)")
    
    # Display results if available
    if hasattr(st.session_state, 'screening_results') and st.session_state.screening_results:
        if not structure_only:
            display_results_table(st.session_state.screening_results)
        with st.container():
            if structure_only:
                display_structure_only_3d_viewer(st.session_state.screening_results, project_name)
                create_visualizations(st.session_state.screening_results, structure_only=True)
            else:
                create_visualizations(st.session_state.screening_results, structure_only=False)
    
    # End of main content and sidebar logic



if __name__ == "__main__":
    main()
