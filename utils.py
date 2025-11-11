from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import os
import subprocess
import yaml
from collections import OrderedDict

def find_similar_molecules(smiles: str, threshold: int):
    columns = ['molecule_chembl_id', 'similarity', 'pref_name', 'molecule_structures']
    try:
        from chembl_webresource_client.new_client import new_client as ch
        res = ch.similarity.filter(smiles=smiles, similarity=threshold).only(columns)
        results = list(res)
        return results
    except Exception as e:
        # Silently handle similarity search errors (e.g., invalid SMILES)
        # Don't display error messages to avoid cluttering the UI
        return None

def get_molecule_by_chembl_id(chembl_id: str):
    """Fetch a molecule by its ChEMBL ID.
    
    Args:
        chembl_id (str): ChEMBL molecule ID (e.g., 'CHEMBL25')
        
    Returns:
        dict: Molecule data with structure and metadata, or None if not found
    """
    try:
        # Validate input
        if not chembl_id or not isinstance(chembl_id, str):
            return None
        
        # Normalize ChEMBL ID (uppercase)
        chembl_id = chembl_id.upper().strip()
        
        # Basic format validation
        if not chembl_id.startswith('CHEMBL'):
            return None
        
        from chembl_webresource_client.new_client import new_client as ch
        
        # Fetch molecule by ChEMBL ID
        molecules = list(ch.molecule.filter(molecule_chembl_id=chembl_id).only([
            'molecule_chembl_id', 'pref_name', 'molecule_structures', 'molecule_properties'
        ]))
        
        if molecules:
            return molecules[0]
        else:
            return None
            
    except Exception as e:
        # Silently handle errors (e.g., invalid ChEMBL ID, network issues)
        return None

def get_smiles_from_molecule(molecule_data):
    """Extract SMILES string from molecule data with fallback options.
    
    Args:
        molecule_data (dict): Molecule data from ChEMBL
        
    Returns:
        str: SMILES string or empty string if not found
    """
    if not molecule_data:
        return ""
    
    structures = molecule_data.get("molecule_structures", {})
    if not structures:
        return ""
    
    # Try different SMILES formats in order of preference
    smiles_formats = [
        "canonical_smiles",  # Most preferred
        "smiles",           # Alternative
        "molfile",          # MOL format (can be converted)
        "inchi"             # InChI format (can be converted)
    ]
    
    for format_name in smiles_formats:
        smiles = structures.get(format_name, "")
        if smiles and isinstance(smiles, str) and smiles.strip():
            return smiles.strip()
    
    return ""

def render_similarity_table(similar_molecules):
    # Example: Render a simple HTML table
    if not similar_molecules:
        return "<p>No similar molecules found.</p>"
    html = "<table border='1'><tr><th>ChEMBL ID</th><th>Similarity</th><th>Name</th><th>Structure</th></tr>"
    for mol in similar_molecules:
        html += f"<tr><td>{mol.get('molecule_chembl_id')}</td><td>{mol.get('similarity')}</td><td>{mol.get('pref_name')}</td><td>{mol.get('molecule_structures', {}).get('canonical_smiles', '')}</td></tr>"
    html += "</table>"
    return html

# Target prediction functions
import onnxruntime
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

FP_SIZE = 1024
RADIUS = 2

# Global variable for the model session
_ort_session = None

def get_ort_session():
    """Get or create the ONNX runtime session for target prediction."""
    global _ort_session
    if _ort_session is None:
        try:
            model_path = "chembl_32_multitask.onnx"
            _ort_session = onnxruntime.InferenceSession(model_path)
        except Exception as e:
            import streamlit as st
            st.warning(f"Could not load target prediction model: {e}")
            return None
    return _ort_session

def calc_morgan_fp(smiles):
    """Calculate Morgan fingerprint for a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
            mol, RADIUS, nBits=FP_SIZE)
        a = np.zeros((0,), dtype=np.float32)
        Chem.DataStructs.ConvertToNumpyArray(fp, a)
        return a
    except Exception as e:
        import streamlit as st
        st.warning(f"Error calculating fingerprint: {e}")
        return None

def format_preds(preds, targets):
    """Format predictions with target names."""
    try:
        preds = np.concatenate(preds).ravel()
        np_preds = [(tar, pre) for tar, pre in zip(targets, preds)]
        dt = [('chembl_id', '|U20'), ('pred', '<f4')]
        np_preds = np.array(np_preds, dtype=dt)
        np_preds[::-1].sort(order='pred')
        return np_preds
    except Exception as e:
        import streamlit as st
        st.warning(f"Error formatting predictions: {e}")
        return None

def predict_targets(smiles):
    """Predict targets for a single SMILES string."""
    try:
        # Get the model session
        ort_session = get_ort_session()
        if ort_session is None:
            return None
            
        # Calculate the fingerprints
        descs = calc_morgan_fp(smiles)
        if descs is None:
            return None

        # Run the prediction
        ort_inputs = {ort_session.get_inputs()[0].name: descs}
        preds = ort_session.run(None, ort_inputs)

        # Format the predictions
        return format_preds(preds, [o.name for o in ort_session.get_outputs()])
    except Exception as e:
        import streamlit as st
        st.warning(f"Error in target prediction: {e}")
        return None

def predict_targets_batch(smiles_list):
    """Predict targets for multiple SMILES strings."""
    preds = []
    for smile in smiles_list:
        pred = predict_targets(smile)
        if pred is not None:
            preds.append(pred)
    if preds:
        return np.concatenate(preds)
    return None

def get_top_targets(predictions, top_n=10, filter_herg=False, admetica_herg=None):
    """Get the top N predicted targets with their scores.
    
    Args:
        predictions: Target predictions array
        top_n: Number of top targets to return
        filter_herg: If True, filter out hERG (CHEMBL240) from results when ADMETICA is available
        admetica_herg: ADMETICA hERG inhibition probability (0-1) to replace target prediction
    """
    if predictions is None:
        return []
    
    try:
        # Convert to list of tuples for easier handling
        target_list = [(str(pred['chembl_id']), float(pred['pred'])) for pred in predictions]
        
        # Handle hERG based on ADMETICA availability
        if admetica_herg is not None:
            # Replace hERG target prediction with ADMETICA prediction
            target_list = [(target_id, score) for target_id, score in target_list if target_id != "CHEMBL240"]
            # Add ADMETICA hERG prediction
            target_list.append(("CHEMBL240", admetica_herg))
        elif filter_herg:
            # Filter out hERG if requested (when ADMETICA is available but no value provided)
            target_list = [(target_id, score) for target_id, score in target_list if target_id != "CHEMBL240"]
        
        # Sort by prediction score (descending)
        target_list.sort(key=lambda x: x[1], reverse=True)
        # Return top N
        return target_list[:top_n]
    except Exception as e:
        import streamlit as st
        st.warning(f"Error getting top targets: {e}")
        return []

def get_target_name(target_id):
    """Get the common name for a ChEMBL target ID."""
    try:
        from chembl_webresource_client.new_client import new_client as ch
        target = ch.target.filter(target_chembl_id=target_id).only(['pref_name', 'organism'])[0]
        if target.get('pref_name'):
            return target['pref_name']
        else:
            return f"Target {target_id}"
    except Exception:
        return f"Target {target_id}"

# Synthetic accessibility and ADMET prediction functions
def calculate_synthetic_accessibility(smiles):
    """Calculate synthetic accessibility score using RDKit's SA_Score."""
    try:
        from rdkit import Chem
        from rdkit.Contrib.SA_Score import sascorer
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Calculate SA score (lower is easier to synthesize)
        sa_score = sascorer.calculateScore(mol)
        return sa_score
    except ImportError:
        import streamlit as st
        st.warning("SA_Score not available. Please install rdkit-contrib: pip install rdkit-contrib")
        return None
    except Exception as e:
        import streamlit as st
        st.warning(f"Error calculating synthetic accessibility: {e}")
        return None

# Boltz integration functions
import json
import re

def parse_protein_chains(protein_sequence):
    """
    Parse protein sequence into individual chains.
    Converts input to uppercase before splitting and storing.
    
    Args:
        protein_sequence (str): Protein sequence(s) separated by colons
        
    Returns:
        list: List of dictionaries with chain information
    """
    protein_sequence = protein_sequence.upper()
    protein_chains = []
    if protein_sequence.strip():
        # Split by colon to get individual chains
        chains = protein_sequence.split(':')
        
        # Validate number of chains (max 23: A-W, X reserved for ligand)
        if len(chains) > 23:
            import streamlit as st
            st.error(f"Too many protein chains ({len(chains)}). Maximum allowed is 23 chains (A-W). Chain X is reserved for ligands.")
            return []
        
        for i, chain in enumerate(chains):
            chain_id = chr(65 + i)  # A, B, C, ..., W (65-87)
            protein_chains.append({
                "protein": {
                    "id": chain_id,
                    "sequence": chain.strip()
                }
            })
    return protein_chains

def create_boltz_yaml(workspace_name, design_name, protein_sequence, ligand_smiles, binding_pocket_constraints=None, cofactor_info=None, template_cif_path=None, ptm_modifications=None):
    """Create a YAML file for Boltz prediction."""
    # Create boltzomics_results directory if it doesn't exist
    results_dir = "boltzomics_results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    # Create filename
    filename = f"{workspace_name}_{design_name}.yaml"
    filepath = os.path.join(results_dir, filename)
    
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
    
    # Create YAML content with protein chains and main ligand
    yaml_content = {
        "sequences": protein_chains + [
            {
                "ligand": {
                    "id": "X",
                    "smiles": ligand_smiles
                }
            }
        ],
    }
    
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
    # Backward compatibility: handle single cofactor dict
    elif cofactor_info and isinstance(cofactor_info, dict) and (cofactor_info.get('smiles') or cofactor_info.get('ccd')):
        cofactor_entry = {
            "ligand": {
                "id": "T"  # Use T for backward compatibility
            }
        }
        
        # Add either SMILES or CCD code
        if cofactor_info.get('smiles'):
            cofactor_entry["ligand"]["smiles"] = cofactor_info['smiles']
        elif cofactor_info.get('ccd'):
            cofactor_entry["ligand"]["ccd"] = cofactor_info['ccd']
        
        # Add co-factor to sequences
        yaml_content["sequences"].append(cofactor_entry)

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
        import copy
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
            # Always append the properties block for affinity prediction
            f.write("properties:\n")
            f.write("  - affinity:\n")
            f.write("      binder: X\n")
    else:
        with open(filepath, 'w') as f:
            yaml.dump(yaml_content, f, default_flow_style=False)
            # Always append the properties block for affinity prediction
            f.write("properties:\n")
            f.write("  - affinity:\n")
            f.write("      binder: X\n")

    # Add templates section if template_cif_path is provided
    if template_cif_path:
        yaml_content["templates"] = [{"cif": template_cif_path}]

    return filepath

def run_boltz_prediction(yaml_filepath, use_gpu=True, override=False, recycling_steps=3, sampling_steps=200, diffusion_samples=1, max_parallel_samples=5, step_scale=1.638, affinity_mw_correction=False, max_msa_seqs=8192, sampling_steps_affinity=200, diffusion_samples_affinity=5, subsample_msa=False, num_subsampled_msa=1024, timeout=300):
    """Run Boltz prediction on the YAML file."""
    try:
        # Change to the directory containing the YAML file
        yaml_dir = os.path.dirname(yaml_filepath)
        yaml_filename = os.path.basename(yaml_filepath)
        # Run boltz predict command
        cmd = ["boltz", "predict", yaml_filename, "--use_msa_server", "--output_format", "pdb"]
        # Add override flag if specified
        if override:
            cmd.append("--override")
        # Add CPU accelerator flag if GPU is disabled
        if not use_gpu:
            cmd.append("--accelerator")
            cmd.append("cpu")
        # Add Boltz parameters
        cmd.extend(["--recycling_steps", str(int(recycling_steps))])
        cmd.extend(["--sampling_steps", str(int(sampling_steps))])
        cmd.extend(["--diffusion_samples", str(int(diffusion_samples))])
        cmd.extend(["--max_parallel_samples", str(int(max_parallel_samples))])
        cmd.extend(["--step_scale", str(float(step_scale))])
        if affinity_mw_correction:
            cmd.append("--affinity_mw_correction")
        cmd.extend(["--max_msa_seqs", str(int(max_msa_seqs))])
        cmd.extend(["--sampling_steps_affinity", str(int(sampling_steps_affinity))])
        cmd.extend(["--diffusion_samples_affinity", str(int(diffusion_samples_affinity))])
        # --- NEW: MSA Subsampling ---
        if subsample_msa:
            cmd.append("--subsample_msa")
            cmd.extend(["--num_subsampled_msa", str(int(num_subsampled_msa))])
        # Print the command line for user reference
        print("[DEBUG]", " ".join(cmd))
        result = subprocess.run(
            cmd,
            cwd=yaml_dir,
            capture_output=True,
            text=True,
            timeout=int(timeout) if timeout else 300
        )
        if result.returncode != 0:
            raise Exception(f"Boltz command failed: {result.stderr}")
        return result.stdout
    except subprocess.TimeoutExpired:
        raise Exception(f"Boltz prediction timed out after {timeout or 300} seconds")
    except Exception as e:
        raise Exception(f"Error running Boltz prediction: {str(e)}")

def parse_boltz_results(yaml_filepath, structure_only=False):
    """Parse Boltz results from the JSON output files."""
    try:
        # Construct the path to the results files
        yaml_dir = os.path.dirname(yaml_filepath)
        yaml_filename = os.path.basename(yaml_filepath)
        yaml_name = os.path.splitext(yaml_filename)[0]
        
        # Path to the affinity results JSON file - corrected path structure with underscores
        affinity_results_path = os.path.join(yaml_dir, f"boltz_results_{yaml_name}", "predictions", yaml_name, f"affinity_{yaml_name}.json")
        
        # Path to the confidence results JSON file - corrected path structure with underscores
        confidence_results_path = os.path.join(yaml_dir, f"boltz_results_{yaml_name}", "predictions", yaml_name, f"confidence_{yaml_name}_model_0.json")
        
        results = {}
        
        if not structure_only:
            # Parse affinity results
            if os.path.exists(affinity_results_path):
                with open(affinity_results_path, 'r') as f:
                    affinity_results = json.load(f)
                # Extract the key values from affinity results
                results.update({
                    "affinity_pred_value": affinity_results.get("affinity_pred_value"),
                    "affinity_probability_binary": affinity_results.get("affinity_probability_binary"),
                    "affinity_pred_value1": affinity_results.get("affinity_pred_value1"),
                    "affinity_probability_binary1": affinity_results.get("affinity_probability_binary1"),
                    "affinity_pred_value2": affinity_results.get("affinity_pred_value2"),
                    "affinity_probability_binary2": affinity_results.get("affinity_probability_binary2")
                })
            else:
                raise Exception(f"Affinity results file not found: {affinity_results_path}")
        # Parse confidence results
        if os.path.exists(confidence_results_path):
            with open(confidence_results_path, 'r') as f:
                confidence_results = json.load(f)
            # Extract the key values from confidence results
            results.update({
                "confidence_score": confidence_results.get("confidence_score"),
                "ptm": confidence_results.get("ptm"),
                "iptm": confidence_results.get("iptm"),
                "ligand_iptm": confidence_results.get("ligand_iptm"),
                "protein_iptm": confidence_results.get("protein_iptm"),
                "complex_plddt": confidence_results.get("complex_plddt"),
                "complex_iplddt": confidence_results.get("complex_iplddt"),
                "complex_pde": confidence_results.get("complex_pde"),
                "complex_ipde": confidence_results.get("complex_ipde"),
                "chains_ptm": confidence_results.get("chains_ptm"),
                "pair_chains_iptm": confidence_results.get("pair_chains_iptm")
            })
        else:
            # If confidence file doesn't exist, use default values
            results.update({
                "confidence_score": 0.0,
                "ptm": 0.0,
                "iptm": 0.0,
                "ligand_iptm": 0.0,
                "protein_iptm": 0.0,
                "complex_plddt": 0.0,
                "complex_iplddt": 0.0,
                "complex_pde": 0.0,
                "complex_ipde": 0.0,
                "chains_ptm": {"0": 0.0, "1": 0.0},
                "pair_chains_iptm": {"0": {"0": 0.0, "1": 0.0}, "1": {"0": 0.0, "1": 0.0}}
            })
        return results
    except Exception as e:
        raise Exception(f"Error parsing Boltz results: {str(e)}")
