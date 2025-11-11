import streamlit as st
import pandas as pd
import numpy as np
import time
import json
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

def create_screening_boltz_yaml(workspace_name, design_name, protein_sequence, ligand_smiles, project_name, binding_pocket_constraints=None, cofactor_info=None, template_cif_path=None, structure_only=False, ptm_modifications=None):
    """Create a YAML file for Boltz prediction in the project-specific directory."""
    # Create project directory
    project_dir = os.path.join("boltzomics_screening_results", project_name)
    if not os.path.exists(project_dir):
        os.makedirs(project_dir)
    
    # Create filename
    filename = f"{workspace_name}_{design_name}.yaml"
    filepath = os.path.join(project_dir, filename)
    
    # Parse protein sequence into chains
    try:
        protein_chains = parse_protein_chains(protein_sequence)
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
    emit_streamlit_feedback=True,
    status_callback: Optional[Callable[[str, Dict[str, Union[str, int, float]]], None]] = None,
):
    """Run Boltz workflow with retry logic and validation."""
    last_error = None
    total_attempts = max_retry_attempts + 1 if enable_retries else 1

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
            )
            utils.run_boltz_prediction(
                yaml_filepath,
                use_gpu,
                override,
                recycling_steps,
                sampling_steps,
                diffusion_samples,
                max_parallel_samples,
                step_scale,
                affinity_mw_correction,
                max_msa_seqs,
                sampling_steps_affinity,
                diffusion_samples_affinity,
                subsample_msa,
                num_subsampled_msa,
                timeout=prediction_timeout_seconds,
            )
            is_valid, validation_error = validate_boltz_results(yaml_filepath, structure_only=structure_only)
            if not is_valid:
                raise Exception(f"Results validation failed: {validation_error}")
            results = utils.parse_boltz_results(yaml_filepath, structure_only=structure_only)
            if results is None:
                raise Exception("Failed to parse Boltz results")
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
    cmd: List[str] = ["boltz", "predict", yaml_filename, "--use_msa_server", "--output_format", "pdb"]
    if params.get("override"):
        cmd.append("--override")
    if not params.get("use_gpu", True):
        cmd.extend(["--accelerator", "cpu"])
    cmd.extend(["--recycling_steps", str(int(params.get("recycling_steps", 3)))])
    cmd.extend(["--sampling_steps", str(int(params.get("sampling_steps", 200)))])
    cmd.extend(["--diffusion_samples", str(int(params.get("diffusion_samples", 1)))])
    cmd.extend(["--max_parallel_samples", str(int(params.get("max_parallel_samples", 5)))])
    cmd.extend(["--step_scale", str(float(params.get("step_scale", 1.638)))])
    if params.get("affinity_mw_correction"):
        cmd.append("--affinity_mw_correction")
    cmd.extend(["--max_msa_seqs", str(int(params.get("max_msa_seqs", 8192)))])
    cmd.extend(["--sampling_steps_affinity", str(int(params.get("sampling_steps_affinity", 200)))])
    cmd.extend(["--diffusion_samples_affinity", str(int(params.get("diffusion_samples_affinity", 5)))])
    if params.get("subsample_msa"):
        cmd.append("--subsample_msa")
        cmd.extend(["--num_subsampled_msa", str(int(params.get("num_subsampled_msa", 1024)))])
    return " ".join(cmd)


def _extract_workspace_design(yaml_name: str) -> Tuple[str, str]:
    if not yaml_name:
        return "screening", "screening"
    return yaml_name, yaml_name


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
    affinity_pred_value = boltz_results.get("affinity_pred_value", 0)
    affinity_probability_binary = boltz_results.get("affinity_probability_binary", 0)
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
            "recycling_steps": params.get("recycling_steps", 3),
            "sampling_steps": params.get("sampling_steps", 200),
            "diffusion_samples": params.get("diffusion_samples", 1),
            "max_parallel_samples": params.get("max_parallel_samples", 5),
            "step_scale": params.get("step_scale", 1.638),
            "affinity_mw_correction": params.get("affinity_mw_correction", False),
            "max_msa_seqs": params.get("max_msa_seqs", 8192),
            "sampling_steps_affinity": params.get("sampling_steps_affinity", 200),
            "diffusion_samples_affinity": params.get("diffusion_samples_affinity", 5),
            "subsample_msa": params.get("subsample_msa", False),
            "num_subsampled_msa": params.get("num_subsampled_msa", 1024),
            "enable_retries": params.get("enable_retries", True),
            "max_retry_attempts": params.get("max_retry_attempts", 2),
            "retry_delay_base": params.get("retry_delay_base", 5),
            "template_cif_path": params.get("template_cif_path"),
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


def execute_screening_job(job: ScreeningJob) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    params = job.parameters or {}
    ligand_smiles = "" if job.structure_only else job.smiles
    ligand_display_name = "" if job.structure_only else job.drug_name
    start_time = time.time()

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
        emit_streamlit_feedback=False,
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
    yaml_filename = f"{job.workspace_name}_{job.design_name}.yaml"
    metadata = {
        "computation_time_seconds": computation_time,
        "timestamp": timestamp,
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


def ensure_job_manager_executor() -> None:
    if not USE_SCREENING_JOB_QUEUE:
        return
    manager = get_job_manager()
    if manager is None:
        return
    session_state = getattr(st, "session_state", None)
    if session_state is not None and session_state.get(JOB_MANAGER_EXECUTOR_STATE_KEY):
        return
    manager.register_executor(execute_screening_job)
    if session_state is not None:
        session_state[JOB_MANAGER_EXECUTOR_STATE_KEY] = True


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

    for protein_name, protein_seq in protein_sequences:
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
            params["binding_pocket_constraints"] = shared_params.get("binding_pocket_constraints")
            params["ptm_modifications"] = shared_params.get("ptm_modifications")

            if use_existing_results:
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
            workspace_name = _generate_workspace_name()
            design_name = _sanitize_design_name(protein_name, drug_name if drug_name else None)
            job_id = f"{project_name}:{uuid.uuid4().hex}"
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

    active_job = summary.get("active_job")
    if active_job:
        ligand_display = active_job.drug_name or ("No ligand" if active_job.structure_only else "")
        st.info(f"Processing: {active_job.protein_name}{' + ' + ligand_display if ligand_display else ''}")

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
    new_computations = 0
    start_time = time.time()

    for protein_name, protein_seq in protein_sequences:
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
                "override": not use_existing_results,
                "recycling_steps": recycling_steps,
                "sampling_steps": sampling_steps,
                "diffusion_samples": diffusion_samples,
                "max_parallel_samples": max_parallel_samples,
                "step_scale": step_scale,
                "affinity_mw_correction": affinity_mw_correction,
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
                "cofactor_info": cofactor_info,
                "ptm_modifications": ptm_modifications,
                "prediction_timeout_seconds": prediction_timeout_seconds,
            }

            workspace_name = _generate_workspace_name()
            design_name = _sanitize_design_name(protein_name, drug_name if drug_name else None)
            boltz_results = None
            ran_prediction = False
            job_elapsed: Optional[float] = None

            if use_existing_results:
                try:
                    existing_yaml_name = find_existing_screening_results(protein_name, drug_name if not structure_only else "", project_name)
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

            if boltz_results is None:
                ran_prediction = True
                new_computations += 1
                job_start_time = time.time()
                ligand_smiles_value = "" if structure_only else drug_smiles_str
                boltz_results, success, error_message = run_boltz_with_retry(
                    workspace_name=workspace_name,
                    design_name=design_name,
                    protein_sequence=protein_seq,
                    ligand_smiles=ligand_smiles_value,
                    project_name=project_name,
                    protein_display_name=protein_name,
                    ligand_display_name=drug_name if not structure_only else "",
                    use_gpu=use_gpu,
                    binding_pocket_constraints=binding_pocket_constraints,
                    override=command_params["override"],
                    recycling_steps=recycling_steps,
                    sampling_steps=sampling_steps,
                    diffusion_samples=diffusion_samples,
                    max_parallel_samples=max_parallel_samples,
                    step_scale=step_scale,
                    affinity_mw_correction=affinity_mw_correction,
                    max_msa_seqs=max_msa_seqs,
                    sampling_steps_affinity=sampling_steps_affinity,
                    diffusion_samples_affinity=diffusion_samples_affinity,
                    cofactor_info=cofactor_info,
                    enable_retries=enable_retries,
                    max_retry_attempts=max_retry_attempts,
                    retry_delay_base=retry_delay_base,
                    subsample_msa=subsample_msa,
                    num_subsampled_msa=num_subsampled_msa,
                    template_cif_path=template_cif_path,
                    structure_only=structure_only,
                    ptm_modifications=ptm_modifications,
                    prediction_timeout_seconds=prediction_timeout_seconds,
                )
                job_elapsed = time.time() - job_start_time
                if not success:
                    st.error(f"All prediction attempts failed for {job_label}: {error_message[:100]}")
                    boltz_results = None

            command_params["cofactor_info"] = cofactor_info
            result_drug_name = "" if structure_only else drug_name
            result_smiles = "" if structure_only else drug_smiles_str

            if boltz_results:
                result = _create_result_entry(
                    protein_name=protein_name,
                    protein_sequence=protein_seq,
                    drug_name=result_drug_name,
                    smiles=result_smiles,
                    workspace_name=workspace_name,
                    design_name=design_name,
                    boltz_results=boltz_results,
                    params=command_params,
                )
            else:
                result = {
                    "protein_name": protein_name,
                    "drug_name": result_drug_name,
                    "protein_sequence": protein_seq,
                    "smiles": result_smiles,
                    "ic50_um": None,
                    "pic50": None,
                    "affinity_probability": None,
                    "confidence": None,
                    "ptm": None,
                    "iptm": None,
                    "avg_plddt": None,
                    "status": "Failed - All retry attempts exhausted",
                    "workspace": workspace_name,
                    "design": design_name,
                    "cofactor_info": cofactor_info,
                    "boltz2_parameters": {
                        "use_gpu": use_gpu,
                        "recycling_steps": recycling_steps,
                        "sampling_steps": sampling_steps,
                        "diffusion_samples": diffusion_samples,
                        "max_parallel_samples": max_parallel_samples,
                        "step_scale": step_scale,
                        "affinity_mw_correction": affinity_mw_correction,
                        "max_msa_seqs": max_msa_seqs,
                        "sampling_steps_affinity": sampling_steps_affinity,
                        "diffusion_samples_affinity": diffusion_samples_affinity,
                        "subsample_msa": subsample_msa,
                        "num_subsampled_msa": num_subsampled_msa,
                        "enable_retries": enable_retries,
                        "max_retry_attempts": max_retry_attempts,
                        "retry_delay_base": retry_delay_base,
                        "template_cif_path": template_cif_path,
                        "override": command_params["override"],
                        "prediction_timeout_seconds": prediction_timeout_seconds,
                    },
                }

            results.append(result)

            if ran_prediction:
                yaml_filename = f"{workspace_name}_{design_name}.yaml"
                boltz_command = _build_boltz_command(yaml_filename, command_params)
                _persist_result_entry(
                    project_name,
                    result,
                    job_elapsed,
                    binding_pocket_constraints=binding_pocket_constraints,
                    template_cif_path=template_cif_path,
                    boltz_command=boltz_command,
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
        "confidence", "ptm", "iptm", "avg_plddt",
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
            "confidence": st.column_config.NumberColumn("Confidence", format="%.3f"),
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
    
    # Sidebar with logo and navigation
    with st.sidebar:
        # st.image(os.path.join("static", "boltzomics", "boltzomics_light.png"), use_container_width=True)

        st.header(":material/settings: Configuration")
        
        # Computation Settings
        st.subheader("Computation Mode")
        use_gpu = st.toggle("Use GPU", value=True, help="Enable GPU acceleration for significantly faster computation. Disable if GPU is unavailable or for CPU-only processing.")
        use_existing_results = st.toggle("Use Existing Results", value=True, help="If enabled, loads previously computed results to save time. If disabled, forces re-computation ensuring fresh results.")
        max_parallel_samples = st.number_input("Max Parallel Samples", min_value=1, value=5, help="Sets the maximum number of samples processed simultaneously. Higher values speed up computation but require more memory.")

        # Sampling Settings
        st.subheader("Structural Sampling")
        recycling_steps = st.number_input("Recycling Steps", min_value=1, value=4, help="Sets the number of iterative refinement steps to improve prediction accuracy. Higher values enhance quality but increase computation time.")
        sampling_steps = st.number_input("Sampling Steps", min_value=1, value=300, help="Defines the number of steps for sampling the model's distribution. More steps improve prediction stability but require more computation.")
        diffusion_samples = st.number_input("Diffusion Samples", min_value=1, value=1, help="Specifies the number of diffusion samples generated per prediction. Increasing this improves robustness but increases runtime.")
        step_scale = st.number_input("Step Scale", value=1.638, format="%.3f", help="Controls the sampling temperature of the diffusion process. Lower values (e.g., 1-1.5) increase diversity, while higher values (e.g., 1.5-2) prioritize precision.")

        # Affinity Prediction Settings
        st.subheader("Affinity Prediction")
        affinity_mw_correction = st.toggle("Molecular Weight Correction", value=False, help="If enabled, applies a molecular weight correction to affinity predictions, improving accuracy for certain molecules. Disable for standard predictions.")
        sampling_steps_affinity = st.number_input("Sampling Steps (Affinity)", min_value=1, value=300, help="Number of sampling steps for affinity predictions. More steps enhance accuracy but extend computation time.")
        diffusion_samples_affinity = st.number_input("Diffusion Samples (Affinity)", min_value=1, value=7, help="Number of diffusion samples for affinity predictions. Higher values improve reliability but increase runtime.")

        # MSA (Multiple Sequence Alignment) Settings
        st.subheader("Multiple Sequence Alignment")
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
                value=st.session_state.get('structure_only', False),
                help="If enabled, skip ligand input and affinity prediction. Only protein structure will be predicted."
            )
        # Store in session state for downstream use
        st.session_state["structure_only"] = structure_only
        
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
                    ensure_job_manager_executor()
                    shared_params = {
                        "use_gpu": use_gpu,
                        "override": not use_existing_results,
                        "recycling_steps": recycling_steps,
                        "sampling_steps": sampling_steps,
                        "diffusion_samples": diffusion_samples,
                        "max_parallel_samples": max_parallel_samples,
                        "step_scale": step_scale,
                        "affinity_mw_correction": affinity_mw_correction,
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
                        "cofactor_info": cofactor_info,
                        "ptm_modifications": ptm_modifications,
                        "prediction_timeout_seconds": prediction_timeout_minutes * 60,
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
                            use_existing_results=use_existing_results,
                            recycling_steps=recycling_steps,
                            sampling_steps=sampling_steps,
                            diffusion_samples=diffusion_samples,
                            max_parallel_samples=max_parallel_samples,
                            step_scale=step_scale,
                            affinity_mw_correction=affinity_mw_correction,
                            max_msa_seqs=max_msa_seqs,
                            sampling_steps_affinity=sampling_steps_affinity,
                            diffusion_samples_affinity=diffusion_samples_affinity,
                            cofactor_info=cofactor_info,
                            binding_pocket_constraints=binding_pocket_constraints,
                            enable_retries=enable_retries,
                            max_retry_attempts=max_retry_attempts,
                            retry_delay_base=retry_delay_base,
                            subsample_msa=subsample_msa,
                            num_subsampled_msa=num_subsampled_msa,
                            template_cif_path=template_cif_path,
                            structure_only=structure_only,
                            ptm_modifications=ptm_modifications,
                            prediction_timeout_seconds=prediction_timeout_minutes * 60,
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
