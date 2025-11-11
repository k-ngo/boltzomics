import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
import math
import os
import zipfile
import io
from datetime import datetime

try:
    from variables import METRIC_DESCRIPTIONS  # type: ignore
    _VARIABLES_IMPORT_ERROR = None
except Exception as _variables_exc:  # pragma: no cover - optional dependency
    METRIC_DESCRIPTIONS = {}
    _VARIABLES_IMPORT_ERROR = str(_variables_exc)

# Universal plot font sizes
UNIVERSAL_AXIS_LABEL_SIZE = 18
UNIVERSAL_TICK_SIZE = 18
UNIVERSAL_TITLE_SIZE = 22


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
    view = mv.view(width=800, height=height, panel=True)
    view.addModel(pdb_data, name=f"{workspace_name}_{design_name}")
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
            selected_pose = None
            for display_name, pose in pose_options:
                if display_name == selected_pose_display:
                    selected_pose = pose
                    break
            
            if selected_pose:
                if not structure_only:
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Protein", selected_pose['protein_name'])
                        st.metric("pIC50", f"{selected_pose['pic50']:.2f}")
                    with col2:
                        st.metric("Drug", selected_pose['drug_name'])
                        st.metric("Confidence", f"{selected_pose['confidence']:.3f}")
                
                
                # Individual PDB download
                pdb_filepath = find_batch_boltz_structure_file(selected_pose['workspace'], selected_pose['design'], project_name)
                if pdb_filepath and os.path.exists(pdb_filepath):
                    file_name = f"{selected_pose['protein_name']}_{selected_pose['drug_name']}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdb"
                    with open(pdb_filepath, "rb") as f:
                        st.download_button(
                            label="Download Current PDB",
                            data=f,
                            file_name=file_name,
                            key=f"download_pdb_{hash(pdb_filepath)}",
                            use_container_width=True
                        )
                else:
                    st.info("PDB file not found for download.")
                
                # Download all PDBs functionality
                if len(available_poses) > 1:
                    # Create zip file in memory
                    zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                        added_files = 0
                        for pose in available_poses:
                            pdb_path = find_batch_boltz_structure_file(pose['workspace'], pose['design'], project_name)
                            if pdb_path and os.path.exists(pdb_path):
                                # Create a descriptive filename
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
                            use_container_width=True
                        )
                    else:
                        st.error("No PDB files found to download.")
        
        with right_col:
            # Display 3D structure
            if selected_pose:
                with st.container(border=False):
                    display_screening_3d_structure(
                        selected_pose['workspace'],
                        selected_pose['design'],
                        project_name,
                        height=600
                    )
            else:
                st.info("Select a pose from the dropdown to view the 3D structure.")
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
                "pdb_filepath": pdb_filepath
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
