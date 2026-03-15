import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
import math
import os
import json
import zipfile
import io
from datetime import datetime
import re

try:
    from variables import METRIC_DESCRIPTIONS  # type: ignore
    _VARIABLES_IMPORT_ERROR = None
except Exception as _variables_exc:  # pragma: no cover - optional dependency
    METRIC_DESCRIPTIONS = {}
    _VARIABLES_IMPORT_ERROR = str(_variables_exc)

try:
    from structure_refinement import (  # type: ignore
        quick_interface_check,
        compare_wt_mutant_quick,
        check_openmm_available,
        minimize_structure_openmm,
    )
except Exception:
    quick_interface_check = None
    compare_wt_mutant_quick = None
    check_openmm_available = None
    minimize_structure_openmm = None

# Universal plot font sizes
UNIVERSAL_AXIS_LABEL_SIZE = 18
UNIVERSAL_TICK_SIZE = 18
UNIVERSAL_TITLE_SIZE = 22
CONFIDENCE_HELP_TEXT = (
    "Confidence = 0.8 x pLDDT + 0.2 x ipTM. "
    "For this formula, pLDDT is normalized to 0-1 (displayed Avg pLDDT is 0-100)."
)


def find_screening_boltz_structure_file(workspace_name, design_name, project_name):
    """Find the Boltz-generated structure file (PDB format) for screening results."""
    # First try with the exact workspace_name from the JSON
    filename = f"{workspace_name}_{design_name}_model_0.pdb"
    boltz_results_dir = f"boltz_results_{workspace_name}_{design_name}"
    predictions_dir = f"{workspace_name}_{design_name}"
    
    possible_paths = [
        # In the correct Boltz output structure for screening (primary location)
        os.path.join("boltzomics_screening_results", project_name, boltz_results_dir, "predictions", predictions_dir, filename),
        # Direct in project directory (fallback)
        os.path.join("boltzomics_screening_results", project_name, filename)
    ]
    
    for filepath in possible_paths:
        if os.path.exists(filepath):
            return filepath
    
    # If exact match fails, try to find files with pattern matching
    # This handles cases where the workspace name in JSON doesn't match actual file names
    project_dir = os.path.join("boltzomics_screening_results", project_name)
    if os.path.exists(project_dir):
        # Look for any boltz_results directory that contains the design_name
        for item in os.listdir(project_dir):
            if item.startswith("boltz_results_") and design_name in item:
                # Extract the actual workspace name from the directory
                # Format: boltz_results_screening_TIMESTAMP_DESIGN_NAME
                try:
                    actual_workspace = item.replace("boltz_results_", "").replace(f"_{design_name}", "")
                    actual_filename = f"{actual_workspace}_{design_name}_model_0.pdb"
                    actual_predictions_dir = f"{actual_workspace}_{design_name}"
                    
                    fallback_path = os.path.join(project_dir, item, "predictions", actual_predictions_dir, actual_filename)
                    if os.path.exists(fallback_path):
                        return fallback_path
                except:
                    continue
    
    return None




def find_batch_boltz_structure_file(workspace_name, design_name, project_name):
    """Backward-compatible wrapper for locating structure files in batch screening outputs."""
    return find_screening_boltz_structure_file(workspace_name, design_name, project_name)


def get_available_pdb_poses(results, project_name):
    """
    Build a list of available protein/ligand poses that have associated PDB files.

    Args:
        results: Iterable of result dictionaries (typically successful screening results)
        project_name: Project directory to search within

    Returns:
        Sorted list of pose dictionaries enriched with file paths and key metrics.
    """
    if not results or not project_name:
        return []

    poses = []
    seen_keys = set()

    for entry in results:
        workspace = entry.get("workspace")
        design = entry.get("design")
        if not workspace or not design:
            continue
        key = (workspace, design)
        if key in seen_keys:
            continue

        pdb_path = find_batch_boltz_structure_file(workspace, design, project_name)
        if not pdb_path or not os.path.exists(pdb_path):
            continue

        seen_keys.add(key)
        raw_drug_name = entry.get("drug_name")
        if raw_drug_name:
            drug_name = raw_drug_name
        else:
            drug_name = "No ligand" if not entry.get("smiles") else "Unknown"

        poses.append(
            {
                "protein_name": entry.get("protein_name") or "Unknown",
                "drug_name": drug_name,
                "workspace": workspace,
                "design": design,
                "pdb_filepath": pdb_path,
                "pic50": entry.get("pic50"),
                "confidence": entry.get("confidence"),
                "affinity_probability": entry.get("affinity_probability"),
                "ic50_um": entry.get("ic50_um"),
                "ptm": entry.get("ptm"),
                "iptm": entry.get("iptm"),
                "avg_plddt": entry.get("avg_plddt"),
                "cofactor_info": entry.get("cofactor_info"),
            }
        )

    poses.sort(
        key=lambda pose: (
            pose.get("pic50") is not None,
            pose.get("pic50") if pose.get("pic50") is not None else float("-inf"),
            pose.get("confidence") if pose.get("confidence") is not None else float("-inf"),
        ),
        reverse=True,
    )
    return poses


def display_screening_3d_structure(
    workspace_name: str,
    design_name: str,
    project_name: str,
    height: int = 600,
    width: int = 600,
) -> None:
    """
    Render a 3D view of a Boltz screening structure using molview.

    Args:
        workspace_name: Workspace portion of the result identifier.
        design_name: Design portion of the result identifier.
        project_name: Active project directory name.
        height: Viewer height in pixels.
        width: Viewer width in pixels.
    """
    pdb_path = find_screening_boltz_structure_file(workspace_name, design_name, project_name)
    if not pdb_path or not os.path.exists(pdb_path):
        st.warning("Structure file not found for 3D visualization.")
        return

    display_screening_3d_structure_from_path(
        pdb_path=pdb_path,
        model_name=f"{workspace_name}_{design_name}",
        height=height,
        width=width,
    )


def _get_relaxed_structure_path(pdb_path: str) -> str:
    """Return expected minimized PDB path for an original PDB file."""
    base, _ = os.path.splitext(pdb_path)
    return f"{base}_minimized.pdb"


def _get_relaxation_metadata_path(relaxed_pdb_path: str) -> str:
    """Return metadata JSON path for a minimized PDB file."""
    base, _ = os.path.splitext(relaxed_pdb_path)
    return f"{base}_meta.json"


def _load_relaxation_metadata(relaxed_pdb_path: str):
    """Load cached relaxation metadata if available."""
    metadata_path = _get_relaxation_metadata_path(relaxed_pdb_path)
    if not os.path.exists(metadata_path):
        return None
    try:
        with open(metadata_path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except Exception:
        return None


def _run_quick_relaxation(
    pdb_path: str,
    max_iterations: int = 150,
    tolerance: float = 20.0,
) -> dict:
    """
    Run a fast OpenMM minimization for interactive UI usage.

    Defaults are tuned for responsiveness (short run) rather than maximum minimization depth.
    """
    if minimize_structure_openmm is None or check_openmm_available is None:
        return {"ok": False, "error": "structure_refinement_unavailable"}
    if not check_openmm_available():
        return {"ok": False, "error": "openmm_not_available"}

    # If input is already a relaxed file, continue minimization in place.
    if pdb_path.endswith("_minimized.pdb"):
        output_path = pdb_path
    else:
        output_path = _get_relaxed_structure_path(pdb_path)
    result = minimize_structure_openmm(
        pdb_path=pdb_path,
        output_path=output_path,
        max_iterations=int(max_iterations),
        tolerance=float(tolerance),
        add_hydrogens=True,
        platform="CPU",
    )

    if not result.minimization_performed or not result.refined_pdb_path:
        return {"ok": False, "error": "minimization_failed"}

    metadata = {
        "performed": bool(result.minimization_performed),
        "refined_pdb": result.refined_pdb_path,
        "initial_energy": result.initial_energy,
        "final_energy": result.final_energy,
        "energy_change": result.energy_change,
        "rmsd": result.rmsd_from_original,
        "clashes_removed": result.clashes_removed,
        "time_seconds": result.refinement_time_seconds,
        "max_iterations": int(max_iterations),
        "tolerance": float(tolerance),
    }
    metadata_path = _get_relaxation_metadata_path(result.refined_pdb_path)
    try:
        with open(metadata_path, "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2)
    except Exception:
        pass

    return {"ok": True, "path": result.refined_pdb_path, "metadata": metadata}


def display_screening_3d_structure_from_path(
    pdb_path: str,
    model_name: str,
    height: int = 600,
    width: int = 600,
    show_panel: bool = True,
) -> None:
    """
    Render a 3D structure from an explicit PDB path using molview.

    Args:
        pdb_path: Path to PDB file.
        model_name: Label used by molview for this model.
        height: Viewer height in pixels.
        width: Viewer width in pixels.
    """
    if not os.path.exists(pdb_path):
        st.warning("Structure file not found for 3D visualization.")
        return

    try:
        import molview as mv  # type: ignore
    except ImportError:
        st.info("3D visualization requires the `molview` package. Install it to enable inline structures.")
        st.code("pip install molview", language="bash")
        return

    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            pdb_data = handle.read()
    except OSError as exc:
        st.error(f"Unable to read structure file: {exc}")
        return

    if not pdb_data.strip():
        st.warning("Structure file is empty. Unable to display 3D structure.")
        return

    # Create molview viewer with dynamic width
    # The width will be overridden by CSS to be responsive
    view = mv.view(width=800, height=height, panel=bool(show_panel))
    view.addModel(pdb_data, name=model_name)
    view.setColorMode('element')
    view.removeSolvent(True)

    # Get the HTML content and wrap it in a responsive container
    html_content = view._repr_html_()

    # Wrap the content in a responsive div with CSS
    responsive_html = f"""
    <style>
        .molview-responsive {{
            width: 100% !important;
            height: {height}px !important;
            display: flex;
            justify-content: center;
        }}
        .molview-responsive > div {{
            width: 100% !important;
            max-width: 100% !important;
        }}
        .molview-responsive canvas {{
            width: 100% !important;
            height: auto !important;
        }}
    </style>
    <div class="molview-responsive">
        {html_content}
    </div>
    """

    components.html(responsive_html, height=int(height) + 50, scrolling=True)


def deduplicate_results(results: list[dict]) -> list[dict]:
    """
    Remove duplicate entries from results, keeping the most complete and recent entries.
    
    Args:
        results: List of prediction result dictionaries
        
    Returns:
        List of deduplicated results
    """
    if not results:
        return results
    
    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(results)
    
    # Create a composite key for identifying duplicates
    df['composite_key'] = df['protein_name'] + '|' + df['drug_name']
    
    # Group by composite key and keep the best entry
    deduplicated_results = []
    
    for key, group in df.groupby('composite_key'):
        if len(group) == 1:
            # No duplicates, keep as is
            result = group.iloc[0].to_dict()
            result.pop('composite_key', None)  # Remove the temporary key
            deduplicated_results.append(result)
        else:
            # Multiple entries for same protein-drug combination
            # Score each entry based on completeness and recency
            best_entry = None
            best_score = -1
            
            for _, row in group.iterrows():
                score = 0
                
                # Completeness score: more non-null values = higher score
                numeric_fields = ['ic50_um', 'pic50', 'affinity_probability', 'confidence', 'ptm', 'iptm', 'avg_plddt']
                for field in numeric_fields:
                    if field in row and pd.notna(row[field]):
                        score += 1
                
                # Status score: Success > other statuses
                if row.get('status') == 'Success':
                    score += 10
                
                # Recency score: prefer entries with workspace/design info (indicating newer runs)
                if row.get('workspace') and row.get('design'):
                    score += 5
                
                # If scores are equal, prefer the first one (original order)
                if score > best_score:
                    best_score = score
                    best_entry = row
            
            if best_entry is not None:
                result = best_entry.to_dict()
                result.pop('composite_key', None)  # Remove the temporary key
                deduplicated_results.append(result)
    
    return deduplicated_results


def _is_wt_name(name: str) -> bool:
    upper = str(name or "").strip().upper()
    return upper in {"WT", "WILDTYPE", "WILD_TYPE"} or upper.endswith("_WT") or upper.endswith("-WT")


def _extract_mutation_positions_from_name(protein_name: str) -> list[int]:
    if not protein_name:
        return []
    matches = re.findall(r"[A-Za-z](\d+)[A-Za-z]", str(protein_name))
    positions: list[int] = []
    for token in matches:
        try:
            pos = int(token)
        except Exception:
            continue
        if pos > 0 and pos not in positions:
            positions.append(pos)
    return sorted(positions)


def _pick_best_wt_pose(poses: list[dict], drug_name: str) -> dict | None:
    candidates = [
        p for p in poses
        if str(p.get("drug_name", "")) == str(drug_name)
        and _is_wt_name(str(p.get("protein_name", "")))
    ]
    if not candidates:
        return None
    candidates.sort(key=lambda item: float(item.get("confidence") or 0.0), reverse=True)
    return candidates[0]


def create_visualizations(results: list[dict], structure_only: bool = False):
    """
    Create visualizations for the screening prediction results.
    
    Args:
        results: List of prediction result dictionaries
        structure_only: If True, skip ligand interaction toggle and metrics in 3D viewer
    """
    if not results:
        return
    
    # Deduplicate results before processing
    deduplicated_results = deduplicate_results(results)
    successful_results = [r for r in deduplicated_results if r["status"] == "Success" and r["pic50"] is not None]
    if not successful_results:
        st.info("No successful results to visualize.")
        return
    df = pd.DataFrame(successful_results)

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

    if not structure_only:
        has_consensus_values = (
            "affinity_pred_value_consensus" in df.columns
            and pd.to_numeric(df["affinity_pred_value_consensus"], errors="coerce").notna().any()
        )
        source_mode = "raw"
        if has_consensus_values:
            affinity_value_source = st.radio(
                "Affinity Value Source for Plots",
                options=[
                    "Raw single-setting output",
                    "Consensus output (multi-sampling aggregate when available)",
                ],
                index=1,
                horizontal=True,
                help=(
                    "Controls pIC50/IC50 and affinity probability values used in plots. "
                    "Consensus uses the configured aggregate (for example full median) when available."
                ),
            )
            source_mode = "consensus" if affinity_value_source.startswith("Consensus") else "raw"
    else:
        source_mode = "raw"

    display_ic50_values = []
    display_pic50_values = []
    display_prob_values = []
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
        if source_mode == "consensus" and consensus_pred is not None:
            selected_pred = consensus_pred
            selected_prob = consensus_prob
        elif source_mode == "raw" and raw_pred is not None:
            selected_pred = raw_pred
            selected_prob = raw_prob
        elif consensus_pred is not None:
            selected_pred = consensus_pred
            selected_prob = consensus_prob
        elif raw_pred is not None:
            selected_pred = raw_pred
            selected_prob = raw_prob

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

    df.loc[:, "ic50_um"] = display_ic50_values
    df.loc[:, "pic50"] = display_pic50_values
    df.loc[:, "affinity_probability"] = display_prob_values
    
    if structure_only:
        # Only show user input summary
        with st.expander("User input summary", expanded=False, icon=":material/list:"):
            # --- Protein FASTA (show all chains as a single entry) ---
            st.markdown("**Protein FASTA**")
            protein_fasta_lines = []
            truncation_found = False
            if not df.empty and "protein_name" in df.columns and "protein_sequence" in df.columns:
                # Convert columns to strings to avoid unhashable type errors
                df_safe = df.copy()
                df_safe["protein_name"] = df_safe["protein_name"].astype(str)
                df_safe["protein_sequence"] = df_safe["protein_sequence"].astype(str)
                
                for pname, seq in df_safe.drop_duplicates(["protein_name"]).loc[:, ["protein_name", "protein_sequence"]].values:
                    protein_fasta_lines.append(f">{pname}")
                    protein_fasta_lines.append(seq)
                    if seq.endswith("..."):
                        truncation_found = True
                st.code("\n".join(protein_fasta_lines), language="text")
                if truncation_found:
                    st.warning("Protein sequences are truncated to the first 50 residues in the results. Full input is not recoverable from results data.")
            else:
                st.info("No protein sequence data available in results.")
            # --- Cofactor Info ---
            st.markdown("**Cofactor info**")
            if "cofactor_info" in df.columns:
                # Convert lists to strings to make them hashable for unique()
                cofactor_series = df["cofactor_info"].dropna()
                if len(cofactor_series) > 0:
                    # Convert each cofactor info to string representation
                    cofactor_strings = []
                    for cinfo in cofactor_series:
                        if isinstance(cinfo, list):
                            cofactor_strings.append(str(cinfo))
                        else:
                            cofactor_strings.append(str(cinfo))
                    unique_cofactors = list(set(cofactor_strings))  # Use set for unique values
                else:
                    unique_cofactors = []
                
                if len(unique_cofactors) == 0:
                    st.info("No cofactor info recorded.")
                else:
                    for cinfo in unique_cofactors:
                        if cinfo == "N/A":
                            st.info("N/A")
                        elif cinfo == "None" or cinfo == "[]" or cinfo == "{}":
                            st.info("None")
                        else:
                            st.write(cinfo)
            else:
                st.info("N/A")
            # --- Template, Binding Pocket, Boltz Command Info ---
            project_data = getattr(st.session_state, 'loaded_project_data', None)
            if project_data:
                template_cif_path = project_data.get("template_cif_path")
                if template_cif_path:
                    st.markdown(f"**Structural template file**")
                    st.code({template_cif_path})
                binding_pocket = project_data.get("binding_pocket_constraints")
                if binding_pocket and isinstance(binding_pocket, dict) and binding_pocket.get('contacts'):
                    contacts = binding_pocket.get('contacts', [])
                    max_distance = binding_pocket.get('max_distance', None)
                    binder = binding_pocket.get('binder', None)
                    st.markdown("**Binding pocket residues**")
                    if contacts:
                        contact_str = ', '.join([f"{c[0]}{c[1]}" for c in contacts])
                        st.write(f"Residues: {contact_str}")
                    if binder:
                        st.write(f"Binder: {binder}")
                    if max_distance is not None:
                        st.write(f"Max distance: {max_distance} Å")
                boltz_cmds = project_data.get("boltz_commands")
                if boltz_cmds:
                    st.markdown("**Boltz predict command(s)**")
                    for i, cmd in enumerate(boltz_cmds):
                        st.code(cmd, language="bash")
        return
    
    # --- Multi-Protein mode: show all visualizations, then user input summary at end ---
    # Filter successful results
    # Filter successful results
    
    # Convert columns to strings to avoid unhashable type errors
    df_safe = df.copy()
    df_safe["drug_name"] = df_safe["drug_name"].astype(str)
    df_safe["protein_name"] = df_safe["protein_name"].astype(str)
    
    # Get unique drugs
    unique_drugs = sorted(df_safe["drug_name"].unique())
    
    # Create tab names
    tab_names = ["Combined Results"] + [f"Drug: {drug}" for drug in unique_drugs]
    
    # Create tabs
    tabs = st.tabs(tab_names)
    # First tab: Combined Results
    with tabs[0]:
        st.subheader(":material/analytics: Combined Analysis")

        # Add multiselects for filtering
        protein_options = df_safe["protein_name"].unique().tolist()
        drug_options = df_safe["drug_name"].unique().tolist()
        selected_proteins = st.multiselect(
            "Select Protein(s)", protein_options, default=protein_options, key="combined_protein_filter"
        )
        selected_drugs = st.multiselect(
            "Select Drug(s)", drug_options, default=drug_options, key="combined_drug_filter"
        )
        # Filter DataFrame based on selections
        filtered_df = df_safe[
            df_safe["protein_name"].isin(selected_proteins) & df_safe["drug_name"].isin(selected_drugs)
        ]
        if filtered_df.empty:
            st.warning("No data to display for the selected mutant(s) and drug(s). Please adjust your selection.")
            return
        # Use filtered_df for all plots in Combined Analysis
        # New layout: left column (3 plots stacked), right column (2 plots stacked)
        left_col, right_col = st.columns(2)
        with left_col:
            # pIC50 by Drug (violin plot)
            num_drugs = len(filtered_df["drug_name"].unique())
            base_height = 400
            min_height_per_drug = 50
            calculated_height = max(base_height, num_drugs * min_height_per_drug)
            max_height = 800
            final_height = min(calculated_height, max_height)
            fig_violin_drug = px.violin(
                filtered_df, 
                y="pic50", 
                x="drug_name",
                title="pIC50 by Drug",
                labels={"pic50": "pIC50", "drug_name": "Drug"},
                color="drug_name",
                color_discrete_sequence=px.colors.qualitative.Dark2,
                box=True  # Show box inside violin
            )
            fig_violin_drug.update_traces(meanline_visible=True, marker=dict(line=dict(color='black', width=1)))
            fig_violin_drug.update_layout(
                height=final_height,
                font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),  
                title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                plot_bgcolor='white',
                paper_bgcolor='white',
                showlegend=False,  # Remove legend since drug names are on x-axis
                margin=dict(l=80, r=80, t=80, b=80)  # Add margins for labels
            )
            fig_violin_drug.update_xaxes(
                tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                tickangle=-45  # Angle labels to prevent overlap
            )
            fig_violin_drug.update_yaxes(
                tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black")
            )
            if calculated_height > max_height:
                st.warning(f"Drug violin plot height was capped due to large number of drugs ({num_drugs}). Consider filtering data for better visualization.")
            st.plotly_chart(fig_violin_drug, use_container_width=True, key="combined_violin_drug")

            # pIC50 by Protein (violin plot)
            num_proteins = len(filtered_df["protein_name"].unique())
            base_height = 400
            min_height_per_protein = 50
            calculated_height = max(base_height, num_proteins * min_height_per_protein)
            max_height = 800
            final_height = min(calculated_height, max_height)
            fig_violin_protein = px.violin(
                filtered_df, 
                y="pic50", 
                x="protein_name",
                title="pIC50 by Protein",
                labels={"pic50": "pIC50", "protein_name": "Protein"},
                color="protein_name",
                color_discrete_sequence=px.colors.qualitative.Dark2,
                box=True  # Show box inside violin
            )
            fig_violin_protein.update_traces(meanline_visible=True, marker=dict(line=dict(color='black', width=1)))
            fig_violin_protein.update_layout(
                height=final_height,
                font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),  
                title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                plot_bgcolor='white',
                paper_bgcolor='white',
                showlegend=False,  # Remove legend since protein names are on x-axis
                margin=dict(l=80, r=80, t=80, b=80)  # Add margins for labels
            )
            fig_violin_protein.update_xaxes(
                tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                tickangle=-45  # Angle labels to prevent overlap
            )
            fig_violin_protein.update_yaxes(
                tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black")
            )
            if calculated_height > max_height:
                st.warning(f"Protein violin plot height was capped due to large number of proteins ({num_proteins}). Consider filtering data for better visualization.")
            st.plotly_chart(fig_violin_protein, use_container_width=True, key="combined_violin_protein")

            # pIC50 vs Affinity Probability scatter plot
            fig_scatter = px.scatter(
                filtered_df, 
                x="pic50", 
                y="affinity_probability",
                title="pIC50 vs Affinity Probability",
                labels={"pic50": "pIC50", "affinity_probability": "Affinity Probability"},
                color="confidence",  # Use confidence for color instead of affinity probability
                color_continuous_scale="viridis",
                hover_data=["protein_name", "drug_name"]
            )
            fig_scatter.update_traces(marker=dict(size=12, line=dict(width=1, color='black')))
            fig_scatter.update_layout(
                height=400,
                font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),  
                title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
            fig_scatter.update_xaxes(tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"))
            fig_scatter.update_yaxes(tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"))
            fig_scatter.update_coloraxes(
                colorbar_tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                colorbar_title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                colorbar=dict(
                    outlinewidth=1, 
                    outlinecolor='black',
                    title="Confidence",
                    thickness=15,
                    lenmode='fraction',
                    len=0.8,
                ),
                cmin=0,  # Set minimum color scale to 0
                cmax=1   # Set maximum color scale to 1
            )
            st.plotly_chart(fig_scatter, use_container_width=True, key="combined_scatter")

        with right_col:
            # pIC50 Heatmap (if possible)
            if len(filtered_df["protein_name"].unique()) > 1 and len(filtered_df["drug_name"].unique()) > 1:
                pivot_df = filtered_df.pivot(index="protein_name", columns="drug_name", values="pic50")
                num_proteins = len(pivot_df.index)
                num_drugs = len(pivot_df.columns)
                base_height = 400
                base_width = 600
                min_height_per_protein = 35
                calculated_height = max(base_height, num_proteins * min_height_per_protein)
                min_width_per_drug = 90
                calculated_width = max(base_width, num_drugs * min_width_per_drug)
                max_height = 1000
                max_width = 1200
                final_height = min(calculated_height, max_height)
                final_width = min(calculated_width, max_width)
                fig_heatmap = px.imshow(
                    pivot_df,
                    title="pIC50 Heatmap",
                    labels=dict(x="Drug", y="Protein", color="pIC50"),
                    color_continuous_scale="thermal",
                    aspect="auto",
                    zmin=pivot_df.min().min(),
                    zmax=pivot_df.max().max()
                )
                fig_heatmap.update_layout(
                    height=final_height,
                    width=final_width,
                    font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),  
                    title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    margin=dict(l=80, r=80, t=80, b=80)
                )
                fig_heatmap.update_xaxes(
                    tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                    title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='black',
                    tickangle=-45,
                    tickmode='array',
                    ticktext=pivot_df.columns.tolist(),
                    tickvals=list(range(len(pivot_df.columns)))
                )
                fig_heatmap.update_yaxes(
                    tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                    title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                    showgrid=True,
                    gridwidth=1,
                    gridcolor='black',
                    tickmode='array',
                    ticktext=pivot_df.index.tolist(),
                    tickvals=list(range(len(pivot_df.index)))
                )
                fig_heatmap.update_coloraxes(
                    colorbar_tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                    colorbar_title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                    colorbar=dict(outlinewidth=1, outlinecolor='black')
                )
                if calculated_height > max_height or calculated_width > max_width:
                    st.warning(f"Heatmap dimensions were capped due to large number of proteins ({num_proteins}) and drugs ({num_drugs}). Consider filtering data for better visualization.")
                st.plotly_chart(fig_heatmap, use_container_width=True, key="combined_heatmap")
            # pIC50 Values (bar plot)
            all_results = filtered_df.sort_values("pic50", ascending=False)
            all_results = all_results.copy()
            all_results['combination'] = all_results['protein_name'] + ' + ' + all_results['drug_name']
            bar_key_right = f"combined_bar_right_{hash(tuple(all_results['combination']))}"
            fig_bar = px.bar(
                all_results,
                x="pic50",
                y="combination",
                title="pIC50 Values",
                labels={"pic50": "pIC50", "combination": ""},
                orientation="h",
                color="drug_name",  # Color by drug name
                color_discrete_sequence=px.colors.qualitative.Set3  # Use discrete colors
            )
            fig_bar.update_layout(
                height=1000,
                font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                plot_bgcolor='white',
                paper_bgcolor='white',
                showlegend=False  # Remove legend
            )
            fig_bar.update_xaxes(tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"))
            fig_bar.update_yaxes(tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"))
            fig_bar.update_traces(marker=dict(line=dict(width=1, color='black')))
            st.plotly_chart(fig_bar, use_container_width=True, key=bar_key_right)
    
    # Individual drug tabs (only if more than 1 drug)
    if len(unique_drugs) > 1:
        for i, drug in enumerate(unique_drugs):
            with tabs[i + 1]:  # +1 because first tab is combined results
                st.subheader(f":material/token: Analysis for {drug}")
                
                # Filter data for this drug using the safe DataFrame
                drug_df = df_safe[df_safe["drug_name"] == drug]
                
                if len(drug_df) == 0:
                    st.info(f"No successful results for {drug}")
                    continue
                
                # Create visualizations for this drug
                col1, col2 = st.columns(2)
                
                with col1:
                    # pIC50 by Protein for this drug (bar plot, sorted by values)
                    drug_df_sorted = drug_df.sort_values('pic50', ascending=True)
                    
                    # Calculate dynamic height based on number of proteins
                    num_proteins = len(drug_df_sorted)
                    base_height = 400
                    min_height_per_protein = 40
                    calculated_height = max(base_height, num_proteins * min_height_per_protein)
                    max_height = 800
                    final_height = min(calculated_height, max_height)
                    
                    fig_drug_protein = px.bar(
                        drug_df_sorted, 
                        x="pic50", 
                        y="protein_name",
                        title=f"pIC50 by Protein - {drug}",
                        labels={"pic50": "pIC50", "protein_name": ""},
                        orientation="h",
                        color="protein_name"  # Color by protein name (different color for each protein)
                    )
                    fig_drug_protein.update_layout(
                        height=final_height,
                        font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),  
                        title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                        plot_bgcolor='white',
                        paper_bgcolor='white',
                        showlegend=False,  # Remove legend since each bar is a different protein
                        margin=dict(l=80, r=80, t=80, b=80)  # Add margins for labels
                    )
                    fig_drug_protein.update_xaxes(
                        tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                        title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black")
                    )
                    fig_drug_protein.update_yaxes(
                        tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                        title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                        tickmode='array',
                        ticktext=drug_df_sorted['protein_name'].tolist(),
                        tickvals=list(range(len(drug_df_sorted)))
                    )
                    # Add black edges around bars
                    fig_drug_protein.update_traces(marker=dict(line=dict(width=1, color='black')))
                    
                    # Show warning if height was capped
                    if calculated_height > max_height:
                        st.warning(f"Chart height was capped due to large number of proteins ({num_proteins}). Consider filtering data for better visualization.")
                    
                    st.plotly_chart(fig_drug_protein, use_container_width=True)
                
                with col2:
                    # pIC50 vs Affinity Probability scatter plot for this drug
                    fig_drug_scatter = px.scatter(
                        drug_df, 
                        x="pic50", 
                        y="affinity_probability",
                        title=f"pIC50 vs Affinity Probability - {drug}",
                        labels={"pic50": "pIC50", "affinity_probability": "Affinity Probability"},
                        color="confidence",  # Use confidence for color instead of affinity probability
                        color_continuous_scale="viridis",
                        hover_data=["protein_name"],
                        size="confidence"  # Use confidence for marker size instead of affinity probability
                    )
                    fig_drug_scatter.update_traces(marker=dict(size=12, line=dict(width=1, color='black')))
                    fig_drug_scatter.update_layout(
                        height=400,
                        font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),  
                        title_font=dict(size=UNIVERSAL_TITLE_SIZE, color="black"),
                        plot_bgcolor='white',
                        paper_bgcolor='white'
                    )
                    fig_drug_scatter.update_xaxes(tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"))
                    fig_drug_scatter.update_yaxes(tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"))
                    fig_drug_scatter.update_coloraxes(
                        colorbar_tickfont=dict(size=UNIVERSAL_TICK_SIZE, color="black"), 
                        colorbar_title_font=dict(size=UNIVERSAL_AXIS_LABEL_SIZE, color="black"),
                        colorbar=dict(
                            outlinewidth=1, 
                            outlinecolor='black',
                            title="Confidence",
                            thickness=15,
                            lenmode='fraction',
                            len=0.8,
                        ),
                        cmin=0,  # Set minimum color scale to 0
                        cmax=1   # Set maximum color scale to 1
                    )
                    # Add black borders around markers
                    st.plotly_chart(fig_drug_scatter, use_container_width=True)
                
                # Summary statistics for this drug
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Predictions", len(drug_df))
                with col2:
                    st.metric("Average pIC50", f"{drug_df['pic50'].mean():.2f}")
                with col3:
                    st.metric("Max pIC50", f"{drug_df['pic50'].max():.2f}")
                with col4:
                    st.metric("Min pIC50", f"{drug_df['pic50'].min():.2f}")

    # Add 3D Structure Visualization Section
    st.subheader(":material/view_in_ar: 3D Structure Visualization")
    
    # Get project name from session state or results
    project_name = st.session_state.get('loaded_project_name', 'Unknown')
    
    # Get available PDB poses
    available_poses = get_available_pdb_poses(successful_results, project_name)
    
    if available_poses:
        # Create two columns for better layout
        left_col, right_col = st.columns([1, 2])
        selected_pose = None
        selected_pdb = None
        displayed_pdb = None
        relaxed_pdb = None
        relaxed_exists = False
        relaxation_meta = None
        viewer_mode = "Original"
        
        with left_col:      
            # Create dropdown options
            pose_options = []
            for pose in available_poses:
                display_name = f"{pose['protein_name']} + {pose['drug_name']} (pIC50: {pose['pic50']:.2f}, Confidence: {pose['confidence']:.3f})"
                pose_options.append((display_name, pose))
            
            # Create dropdown
            selected_pose_display = st.selectbox(
                "Select a protein-drug combination to visualize its 3D docked structure:" if not structure_only else "Select a protein to visualize its 3D structure:",
                options=[opt[0] for opt in pose_options],
                key="pose_selector" if not structure_only else "structure_only_pose_selector"
            )
            
            # Find the selected pose
            for display_name, pose in pose_options:
                if display_name == selected_pose_display:
                    selected_pose = pose
                    break
            
            if selected_pose:
                if not structure_only:
                    st.markdown(
                        (
                            "<p style='font-size:1.05rem; font-weight:600; margin:0 0 0.4rem 0;'>"
                            f"Protein: {selected_pose['protein_name']} | Drug: {selected_pose['drug_name']}"
                            "</p>"
                        ),
                        unsafe_allow_html=True,
                    )
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        if selected_pose.get("pic50") is not None:
                            st.metric(
                                "pIC50",
                                f"{float(selected_pose['pic50']):.2f}",
                                help=METRIC_DESCRIPTIONS.get("pic50", "Predicted binding potency (higher is better)."),
                            )
                    with col2:
                        if selected_pose.get("confidence") is not None:
                            st.metric(
                                "Confidence",
                                f"{float(selected_pose['confidence']):.3f}",
                                help=CONFIDENCE_HELP_TEXT,
                            )
                    with col3:
                        if selected_pose.get("avg_plddt") is not None:
                            st.metric(
                                "Avg pLDDT",
                                f"{float(selected_pose['avg_plddt']):.1f}",
                                help=METRIC_DESCRIPTIONS.get("plddt", "Average per-residue confidence (0-100)."),
                            )
                # Individual PDB download
                pdb_filepath = find_batch_boltz_structure_file(selected_pose['workspace'], selected_pose['design'], project_name)
                if pdb_filepath and os.path.exists(pdb_filepath):
                    selected_pdb = pdb_filepath
                    relaxed_pdb = _get_relaxed_structure_path(selected_pdb)
                    relaxed_exists = os.path.exists(relaxed_pdb)

                    file_name_orig = (
                        f"{selected_pose['protein_name']}_{selected_pose['drug_name']}_original_"
                        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdb"
                    )
                    with open(selected_pdb, "rb") as f:
                        st.download_button(
                            label="Download Original PDB",
                            data=f.read(),
                            file_name=file_name_orig,
                            key=f"download_original_pdb_{hash(selected_pdb)}",
                            use_container_width=True,
                        )
                    if relaxed_exists:
                        file_name_relaxed = (
                            f"{selected_pose['protein_name']}_{selected_pose['drug_name']}_relaxed_"
                            f"{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdb"
                        )
                        with open(relaxed_pdb, "rb") as f:
                            st.download_button(
                                label="Download Relaxed PDB",
                                data=f.read(),
                                file_name=file_name_relaxed,
                                key=f"download_relaxed_pdb_{hash(relaxed_pdb)}",
                                use_container_width=True,
                            )
                else:
                    st.info("PDB file not found for download.")

                # Optional on-demand post-prediction relaxation.
                selected_pdb = pdb_filepath if pdb_filepath and os.path.exists(pdb_filepath) else None
                if selected_pdb:
                    relaxed_pdb = _get_relaxed_structure_path(selected_pdb)
                    relaxed_exists = os.path.exists(relaxed_pdb)
                    relaxation_meta = _load_relaxation_metadata(relaxed_pdb) if relaxed_exists else None

                    with st.expander("Post-prediction Relaxation.", expanded=False):
                        if check_openmm_available is None:
                            st.caption("OpenMM relaxation is unavailable in this environment.")
                        elif not check_openmm_available():
                            st.caption("OpenMM is not installed. Install it to enable relaxation.")
                        else:
                            c1, c2 = st.columns(2)
                            with c1:
                                relax_iterations = st.slider(
                                    "Quick Relax Iterations",
                                    min_value=50,
                                    max_value=500,
                                    value=150,
                                    step=25,
                                    key=f"relax_iter_{selected_pose['workspace']}_{selected_pose['design']}",
                                    help="Lower values complete faster. 150 is a quick interactive default.",
                                )
                            with c2:
                                relax_tolerance = st.slider(
                                    "Tolerance (kJ/mol/nm)",
                                    min_value=5.0,
                                    max_value=50.0,
                                    value=20.0,
                                    step=1.0,
                                    key=f"relax_tol_{selected_pose['workspace']}_{selected_pose['design']}",
                                    help="Higher tolerance converges faster with less strict minimization.",
                                )

                            relax_input_pdb = relaxed_pdb if relaxed_exists else selected_pdb
                            if st.button(
                                "Relax Model",
                                key=f"run_relax_{selected_pose['workspace']}_{selected_pose['design']}",
                                use_container_width=True,
                                type="tertiary",
                            ):
                                with st.spinner("Running relaxation..."):
                                    relax_result = _run_quick_relaxation(
                                        relax_input_pdb,
                                        max_iterations=relax_iterations,
                                        tolerance=relax_tolerance,
                                    )
                                if relax_result.get("ok"):
                                    if relaxed_exists:
                                        st.success("Relaxation updated.")
                                    else:
                                        st.success("Relaxed structure generated.")
                                    relaxed_pdb = str(relax_result.get("path"))
                                    relaxed_exists = os.path.exists(relaxed_pdb)
                                    relaxation_meta = relax_result.get("metadata")
                                else:
                                    st.error(f"Relaxation failed: {relax_result.get('error', 'unknown_error')}")

                            if relaxed_exists:
                                st.caption("Relaxed structure is available for viewing.")
                                if relaxation_meta:
                                    m1, m2 = st.columns(2)
                                    with m1:
                                        rmsd_val = relaxation_meta.get("rmsd")
                                        st.metric("RMSD", "N/A" if rmsd_val is None else f"{float(rmsd_val):.3f} A")
                                    with m2:
                                        clashes_removed = relaxation_meta.get("clashes_removed")
                                        st.metric("Clashes Removed", int(clashes_removed) if clashes_removed is not None else 0)
                                    runtime_sec = relaxation_meta.get("time_seconds")
                                    runtime_text = (
                                        f"Last relaxation runtime: {float(runtime_sec):.1f}s"
                                        if runtime_sec is not None
                                        else "Last relaxation runtime: N/A"
                                    )
                                    st.caption(runtime_text)
                                if st.button(
                                    "Delete Relaxed Structure",
                                    key=f"delete_relaxed_{selected_pose['workspace']}_{selected_pose['design']}",
                                    use_container_width=True,
                                    type="tertiary",
                                ):
                                    try:
                                        os.remove(relaxed_pdb)
                                    except Exception:
                                        pass
                                    try:
                                        os.remove(_get_relaxation_metadata_path(relaxed_pdb))
                                    except Exception:
                                        pass
                                    st.success("Relaxed structure deleted.")
                                    st.rerun()
                            else:
                                st.caption("No relaxed structure generated yet.")
        
        with right_col:
            # Display 3D structure
            if selected_pose:
                mode_options = ["Original"]
                if relaxed_exists:
                    mode_options.extend(["Relaxed", "Compare"])
                mode_key = f"viewer_mode_{selected_pose['workspace']}_{selected_pose['design']}"
                viewer_mode = st.session_state.get(mode_key, "Original")
                if viewer_mode not in mode_options:
                    viewer_mode = "Original"
                displayed_pdb = relaxed_pdb if viewer_mode == "Relaxed" and relaxed_exists else selected_pdb

                with st.container(border=False):
                    if not selected_pdb or not os.path.exists(selected_pdb):
                        st.info("Selected pose has no available PDB for visualization.")
                    elif viewer_mode == "Compare" and relaxed_pdb and os.path.exists(relaxed_pdb):
                        c1, c2 = st.columns(2)
                        with c1:
                            st.markdown(
                                "<p style='text-align:center; margin:0 0 0.25rem 0; font-weight:600;'>Original</p>",
                                unsafe_allow_html=True,
                            )
                            display_screening_3d_structure_from_path(
                                pdb_path=selected_pdb,
                                model_name=f"{selected_pose['workspace']}_{selected_pose['design']}_original",
                                height=560,
                                show_panel=False,
                            )
                        with c2:
                            st.markdown(
                                "<p style='text-align:center; margin:0 0 0.25rem 0; font-weight:600;'>Relaxed</p>",
                                unsafe_allow_html=True,
                            )
                            display_screening_3d_structure_from_path(
                                pdb_path=relaxed_pdb,
                                model_name=f"{selected_pose['workspace']}_{selected_pose['design']}_relaxed",
                                height=560,
                                show_panel=False,
                            )
                    else:
                        display_screening_3d_structure_from_path(
                            pdb_path=displayed_pdb or selected_pdb,
                            model_name=f"{selected_pose['workspace']}_{selected_pose['design']}",
                            height=600,
                        )

                controls_left, controls_right = st.columns([2, 1])
                with controls_left:
                    st.radio(
                        "Structure Version",
                        options=mode_options,
                        horizontal=True,
                        key=mode_key,
                    )
                with controls_right:
                    if len(available_poses) > 1:
                        zip_buffer = io.BytesIO()
                        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                            added_files = 0
                            for pose in available_poses:
                                pdb_path = find_batch_boltz_structure_file(pose['workspace'], pose['design'], project_name)
                                if pdb_path and os.path.exists(pdb_path):
                                    safe_protein = "".join(c for c in pose['protein_name'] if c.isalnum() or c in (' ', '-', '_')).rstrip()
                                    safe_drug = "".join(c for c in pose['drug_name'] if c.isalnum() or c in (' ', '-', '_')).rstrip()
                                    zip_filename = f"{safe_protein}_{safe_drug}_pIC50_{pose['pic50']:.2f}.pdb"
                                    zip_file.write(pdb_path, zip_filename)
                                    added_files += 1

                        if added_files > 0:
                            zip_buffer.seek(0)
                            zip_filename = f"{project_name}_all_pdbs_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
                            st.download_button(
                                label="Download All PDBs",
                                data=zip_buffer.getvalue(),
                                file_name=zip_filename,
                                mime="application/zip",
                                key="download_all_pdbs_zip",
                                use_container_width=True,
                            )
            else:
                st.info("Select a pose from the dropdown to view the 3D structure.")

        # Dedicated interaction report panel under 3D viewer.
        if selected_pose:
            st.subheader(":material/hub: Protein-Ligand Interaction Report")
            interaction_pdb = displayed_pdb or selected_pdb or find_batch_boltz_structure_file(
                selected_pose.get("workspace"),
                selected_pose.get("design"),
                project_name,
            )
            if not interaction_pdb or not os.path.exists(interaction_pdb):
                st.warning("Structure file not found for interaction analysis.")
            elif quick_interface_check is None:
                st.info(
                    "Interaction analysis backend is unavailable in this environment "
                    "(missing `structure_refinement` dependencies)."
                )
            else:
                try:
                    base_interface = quick_interface_check(interaction_pdb, auto_detect_chains=True)
                except Exception as exc:
                    base_interface = {"error": str(exc)}

                if "error" in base_interface:
                    st.warning(f"Could not analyze selected complex interactions: {base_interface['error']}")
                else:
                    m1, m2, m3, m4, m5 = st.columns(5)
                    with m1:
                        st.metric("Contacts", int(base_interface.get("contacts", 0)))
                        st.caption("Protein-ligand residue contacts within interface cutoff.")
                    with m2:
                        st.metric("H-Bonds", int(base_interface.get("h_bonds", 0)))
                        st.caption("Hydrogen-bond count from interface contact detection.")
                    with m3:
                        st.metric("Hydrophobic", int(base_interface.get("hydrophobic", 0)))
                        st.caption("Hydrophobic contact count at the interface.")
                    with m4:
                        st.metric("Clashes", int(base_interface.get("clashes", 0)))
                        st.caption("Steric overlaps detected at the interface (lower is better).")
                    with m5:
                        clash_score = base_interface.get("clash_score", None)
                        st.metric(
                            "Clash Score",
                            "N/A" if clash_score is None else f"{float(clash_score):.3f}",
                            help="Normalized steric-overlap burden at the interface. Lower is better; values near 0 indicate cleaner geometry.",
                        )
                        st.caption("Normalized steric-overlap burden at the interface. Lower is better.")

                    # Contribution tables for transparency (reviewer-safe, quantitative).
                    with st.expander("Contacts contributors", expanded=False):
                        rows = base_interface.get("contacts_rows", []) or []
                        if rows:
                            contact_df = pd.DataFrame(rows)
                            if "distance_A" in contact_df.columns:
                                contact_df["distance_A"] = pd.to_numeric(contact_df["distance_A"], errors="coerce").round(3)
                            st.dataframe(contact_df, use_container_width=True, hide_index=True)
                        else:
                            st.caption("No contact contributors detected.")

                    with st.expander("H-Bonds contributors", expanded=False):
                        rows = base_interface.get("hbond_rows", []) or []
                        if rows:
                            hbond_df = pd.DataFrame(rows)
                            if "distance_A" in hbond_df.columns:
                                hbond_df["distance_A"] = pd.to_numeric(hbond_df["distance_A"], errors="coerce").round(3)
                            st.dataframe(hbond_df, use_container_width=True, hide_index=True)
                        else:
                            st.caption("No hydrogen-bond contributors detected.")

                    with st.expander("Hydrophobic contributors", expanded=False):
                        rows = base_interface.get("hydrophobic_rows", []) or []
                        if rows:
                            hydro_df = pd.DataFrame(rows)
                            if "distance_A" in hydro_df.columns:
                                hydro_df["distance_A"] = pd.to_numeric(hydro_df["distance_A"], errors="coerce").round(3)
                            st.dataframe(hydro_df, use_container_width=True, hide_index=True)
                        else:
                            st.caption("No hydrophobic contributors detected.")

                    with st.expander("Clash contributors", expanded=False):
                        rows = base_interface.get("clash_rows", []) or []
                        if rows:
                            clash_df = pd.DataFrame(rows)
                            if "distance_A" in clash_df.columns:
                                clash_df["distance_A"] = pd.to_numeric(clash_df["distance_A"], errors="coerce").round(3)
                            st.dataframe(clash_df, use_container_width=True, hide_index=True)
                        else:
                            st.caption("No clash contributors detected.")

                # WT vs mutant interaction comparison for selected drug.
                if compare_wt_mutant_quick is None:
                    st.caption("WT-vs-mutant interaction comparison unavailable (missing dependency).")
                else:
                    selected_protein = str(selected_pose.get("protein_name", ""))
                    selected_drug = str(selected_pose.get("drug_name", ""))
                    mutation_positions = _extract_mutation_positions_from_name(selected_protein)
                    wt_pose = _pick_best_wt_pose(available_poses, selected_drug)

                    if not mutation_positions:
                        st.caption("No mutation code detected in selected protein name; skipping WT-vs-mutant comparison.")
                    elif _is_wt_name(selected_protein):
                        st.caption("Selected pose is WT; choose a mutant pose to see WT-vs-mutant interaction deltas.")
                    elif wt_pose is None:
                        st.caption(
                            f"No WT reference found for drug `{selected_drug}` in current results; "
                            "WT-vs-mutant interaction delta is unavailable."
                        )
                    else:
                        wt_pdb = find_batch_boltz_structure_file(
                            wt_pose.get("workspace"),
                            wt_pose.get("design"),
                            project_name,
                        )
                        if not wt_pdb or not os.path.exists(wt_pdb):
                            st.caption("WT structure file not found for interaction comparison.")
                        else:
                            total_new_clashes = 0
                            disruptions = []
                            similarities = []
                            for pos in mutation_positions:
                                try:
                                    comp = compare_wt_mutant_quick(
                                        wt_pdb=wt_pdb,
                                        mutant_pdb=interaction_pdb,
                                        mutation_position=int(pos),
                                        auto_detect_chains=True,
                                    )
                                except Exception as exc:
                                    comp = {"error": str(exc)}
                                if "error" in comp:
                                    continue
                                disruptions.append(float(comp.get("disruption_score", 0.0)))
                                similarities.append(float(comp.get("fingerprint_similarity", 0.0)))
                                total_new_clashes += int(comp.get("new_clashes", 0))
                            if disruptions:
                                st.markdown("**WT vs Mutant Interaction Delta (same drug)**")
                                d1, d2, d3 = st.columns(3)
                                with d1:
                                    st.metric(
                                        "Mean Disruption",
                                        f"{(sum(disruptions) / len(disruptions)):.3f}",
                                        help="Computed as the mean of per-mutation disruption scores. "
                                        "Per-mutation score = 0.3*(contact-change fraction) + 0.4*(new-clash penalty, capped at 3 clashes) "
                                        "+ 0.3*(1 - fingerprint similarity), with 1.5x boost if mutation is at interface. "
                                        "Lower is better; near 0 means WT-like interface.",
                                    )
                                    st.caption("Weighted contact-change + clash + fingerprint shift score (0-1). Lower is better.")
                                with d2:
                                    st.metric(
                                        "Mean Fingerprint Similarity",
                                        f"{(sum(similarities) / len(similarities)):.3f}",
                                        help="Mean Tanimoto similarity between WT and mutant interaction fingerprints across mutation positions. "
                                        "Fingerprint encodes residue-contact and interaction-type bits. "
                                        "Higher is better; 1.0 means identical fingerprint.",
                                    )
                                    st.caption("Mean Tanimoto similarity of WT-vs-mutant interaction fingerprints. Higher is better.")
                                with d3:
                                    st.metric(
                                        "Total New Clashes",
                                        int(total_new_clashes),
                                        help="Sum of clash events present in mutant but absent in WT across evaluated mutation positions. "
                                        "Lower is better; 0 is ideal.",
                                    )
                                    st.caption("New steric overlaps introduced in mutant vs WT; lower is better.")
                            else:
                                st.caption("WT-vs-mutant interaction comparison could not be computed for detected mutation positions.")
    else:
        st.info("No 3D structures available for visualization. Run Boltz predictions to generate PDB files.")

    # At the end of create_visualizations, after all visualizations and before function end
    # Add variable explanations expander if metadata is available
    if _VARIABLES_IMPORT_ERROR:
        st.info(f"Variable descriptions unavailable ({_VARIABLES_IMPORT_ERROR}).")
    else:
        with st.expander("Explanation of variables", expanded=False, icon=":material/info:"):
            variable_keys = [
                ("pic50", "pIC50"),
                ("predicted_ic50", "IC50 (μM)"),
                ("affinity", "Affinity Probability"),
                ("conf", "Confidence"),
                ("ptm", "pTM"),
                ("iptm", "ipTM"),
                ("plddt", "Avg pLDDT")
            ]
            for key, label in variable_keys:
                desc = METRIC_DESCRIPTIONS.get(key, "No description available.")
                st.markdown(f"- **{label}:** {desc}")
    # Add user input summary expander for reproducibility
    with st.expander("User input summary", expanded=False, icon=":material/list:"):
        # --- Protein FASTA (show all chains as a single entry) ---
        st.markdown("**Protein FASTA**")
        protein_fasta_lines = []
        truncation_found = False
        if not df.empty and "protein_name" in df.columns and "protein_sequence" in df.columns:
            # Convert columns to strings to avoid unhashable type errors
            df_safe = df.copy()
            df_safe["protein_name"] = df_safe["protein_name"].astype(str)
            df_safe["protein_sequence"] = df_safe["protein_sequence"].astype(str)
            
            for pname, seq in df_safe.drop_duplicates(["protein_name"]).loc[:, ["protein_name", "protein_sequence"]].values:
                protein_fasta_lines.append(f">{pname}")
                protein_fasta_lines.append(seq)
                if seq.endswith("..."):
                    truncation_found = True
            st.code("\n".join(protein_fasta_lines), language="text")
            if truncation_found:
                st.warning("Protein sequences are truncated to the first 50 residues in the results. Full input is not recoverable from results data.")
        else:
            st.info("No protein sequence data available in results.")

        # --- Drug SMILES FASTA (full, only if not structure_only) ---
        if not structure_only:
            st.markdown("**Drug SMILES FASTA**")
            smiles_fasta_lines = []
            if not df.empty and "drug_name" in df.columns and "smiles" in df.columns:
                # Convert columns to strings to avoid unhashable type errors
                df_safe = df.copy()
                df_safe["drug_name"] = df_safe["drug_name"].astype(str)
                df_safe["smiles"] = df_safe["smiles"].astype(str)
                
                for dname, smiles in df_safe.drop_duplicates(["drug_name"]).loc[:, ["drug_name", "smiles"]].values:
                    smiles_fasta_lines.append(f">{dname}")
                    smiles_fasta_lines.append(smiles)
                st.code("\n".join(smiles_fasta_lines), language="text")
            else:
                st.info("No drug SMILES data available in results.")

        # --- Cofactor Info ---
        st.markdown("**Cofactor info**")
        if "cofactor_info" in df.columns:
            # Convert lists to strings to make them hashable for unique()
            cofactor_series = df["cofactor_info"].dropna()
            if len(cofactor_series) > 0:
                # Convert each cofactor info to string representation
                cofactor_strings = []
                for cinfo in cofactor_series:
                    if isinstance(cinfo, list):
                        cofactor_strings.append(str(cinfo))
                    else:
                        cofactor_strings.append(str(cinfo))
                unique_cofactors = list(set(cofactor_strings))  # Use set for unique values
            else:
                unique_cofactors = []
            
            if len(unique_cofactors) == 0:
                st.info("No cofactor info recorded.")
            else:
                for cinfo in unique_cofactors:
                    if cinfo == "N/A":
                        st.info("N/A")
                    elif cinfo == "None" or cinfo == "[]" or cinfo == "{}":
                        st.info("None")
                    else:
                        st.write(cinfo)
        else:
            st.info("N/A")

        # --- Template, Binding Pocket, Boltz Command Info ---
        project_data = getattr(st.session_state, 'loaded_project_data', None)
        if project_data:
            # Template file
            template_cif_path = project_data.get("template_cif_path")
            if template_cif_path:
                st.markdown(f"**Structural template file**")
                st.code({template_cif_path})
            # Binding pocket
            binding_pocket = project_data.get("binding_pocket_constraints")
            if binding_pocket and isinstance(binding_pocket, dict) and binding_pocket.get('contacts'):
                contacts = binding_pocket.get('contacts', [])
                max_distance = binding_pocket.get('max_distance', None)
                binder = binding_pocket.get('binder', None)
                st.markdown("**Binding pocket residues**")
                if contacts:
                    contact_str = ', '.join([f"{c[0]}{c[1]}" for c in contacts])
                    st.write(f"Residues: {contact_str}")
                if binder:
                    st.write(f"Binder: {binder}")
                if max_distance is not None:
                    st.write(f"Max distance: {max_distance} Å")
            # Boltz command(s)
            boltz_cmds = project_data.get("boltz_commands")
            if boltz_cmds:
                st.markdown("**Boltz predict command(s)**")
                for i, cmd in enumerate(boltz_cmds):
                    st.code(cmd, language="bash")


def display_structure_only_3d_viewer(results, project_name):
    """
    Display 3D structure visualization for structure-only screening results.
    Args:
        results: List of prediction result dictionaries
        project_name: Name of the project
    """
    if not results:
        st.info("No results to display.")
        return
    # Filter successful results (ignore affinity fields)
    successful_results = [r for r in results if r.get("status") == "Success"]
    if not successful_results:
        st.info("No successful structure predictions to visualize.")
        return
    st.subheader(":material/view_in_ar: 3D Structure Visualization")
    # List available poses (proteins)
    available_poses = []
    for result in successful_results:
        workspace_name = result.get("workspace")
        design_name = result.get("design")
        protein_name = result.get("protein_name", "Unknown")
        # Skip if workspace, design, or project_name is missing
        if not workspace_name or not design_name or not project_name:
            continue
        pdb_filepath = find_batch_boltz_structure_file(workspace_name, design_name, project_name)
        if pdb_filepath:
            available_poses.append({
                "protein_name": protein_name,
                "workspace": workspace_name,
                "design": design_name,
                "pdb_filepath": pdb_filepath,
                "confidence": result.get("confidence"),
                "ptm": result.get("ptm"),
                "iptm": result.get("iptm"),
                "avg_plddt": result.get("avg_plddt"),
            })
    if not available_poses:
        st.info("No 3D structures available for visualization. Run Boltz predictions to generate PDB files.")
        return
    # Dropdown for protein selection
    pose_options = [f"{pose['protein_name']}" for pose in available_poses]
    selected_pose_display = st.selectbox(
        "Select a protein to visualize its 3D structure:",
        options=pose_options,
        key="structure_only_pose_selector"
    )
    selected_pose = next((pose for pose in available_poses if pose['protein_name'] == selected_pose_display), None)
    if selected_pose:
        metric_col1, metric_col2, metric_col3 = st.columns(3)
        with metric_col1:
            if selected_pose.get("confidence") is not None:
                st.metric(
                    "Confidence",
                    f"{float(selected_pose['confidence']):.3f}",
                    help=CONFIDENCE_HELP_TEXT,
                )
        with metric_col2:
            if selected_pose.get("iptm") is not None:
                st.metric(
                    "ipTM",
                    f"{float(selected_pose['iptm']):.3f}",
                    help=METRIC_DESCRIPTIONS.get("iptm", "Predicted interface quality (0-1)."),
                )
        with metric_col3:
            if selected_pose.get("avg_plddt") is not None:
                st.metric(
                    "Avg pLDDT",
                    f"{float(selected_pose['avg_plddt']):.1f}",
                    help=METRIC_DESCRIPTIONS.get("plddt", "Average per-residue confidence (0-100)."),
                )
        # Download button
        pdb_filepath = selected_pose['pdb_filepath']
        if pdb_filepath and os.path.exists(pdb_filepath):
            file_name = f"{selected_pose['protein_name']}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdb"
            with open(pdb_filepath, "rb") as f:
                st.download_button(
                    label="Download Current PDB",
                    data=f,
                    file_name=file_name,
                    key=f"download_pdb_{hash(pdb_filepath)}",
                    icon=":material/download:",
                    use_container_width=True
                )
        # 3D structure viewer
        with st.container(border=False):
            display_screening_3d_structure(
                selected_pose['workspace'],
                selected_pose['design'],
                project_name,
                height=900
            )
    # Download all PDBs as zip
    if len(available_poses) > 1:
        import io, zipfile
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            added_files = 0
            for pose in available_poses:
                pdb_path = pose['pdb_filepath']
                if pdb_path and os.path.exists(pdb_path):
                    safe_protein = "".join(c for c in pose['protein_name'] if c.isalnum() or c in (' ', '-', '_')).rstrip()
                    zip_filename = f"{safe_protein}.pdb"
                    zip_file.write(pdb_path, zip_filename)
                    added_files += 1
        if added_files > 0:
            zip_buffer.seek(0)
            zip_filename = f"{project_name}_all_pdbs_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
            st.download_button(
                label="Download All PDBs",
                data=zip_buffer.getvalue(),
                file_name=zip_filename,
                mime="application/zip",
                key="download_all_pdbs_zip_structure_only",
                icon=":material/archive:",
                use_container_width=True
            )
        else:
            st.error("No PDB files found to download.")
