from typing import List, Tuple, Dict, Optional
import streamlit as st
import re
import pandas as pd
import requests
import time

from mutation_discovery import discover_mutations_for_sequence

@st.dialog("Sequence Viewer - Hover over residues to see indices", width="large")
def show_sequence_viewer_dialog(protein_sequences: List[Tuple[str, str]]):
    """Dialog to show protein sequences with hover tooltips for residue 1-based indices."""
    st.markdown("**Hover over residues to see their 1-based indices**")
    
    if protein_sequences:
        for protein_name, protein_seq in protein_sequences:
            with st.expander(f"**{protein_name}**", expanded=False):
                # Parse chains if multi-chain
                is_valid, _, chains_dict, _ = validate_protein_sequence(protein_seq)
                if is_valid and chains_dict:
                    # Multi-chain protein
                    for chain_id in sorted(chains_dict.keys()):
                        chain_seq = chains_dict[chain_id]
                        st.markdown(f"**Chain {chain_id}:**")
                        
                        # Create sequence display with hover tooltips
                        seq_html = "<div style='font-family: monospace; font-size: 20px; line-height: 1.2;'>"
                        for i, residue in enumerate(chain_seq):
                            # Calculate position in chain (1-based)
                            pos_in_chain = i + 1
                            # Add hover tooltip with residue code and index
                            seq_html += f'<span title="[ {residue} ] Index: {pos_in_chain}" style="cursor: pointer; padding: 1px; border-radius: 2px; font-size: 20px;" onmouseover="this.style.backgroundColor=\'yellow\'; this.style.fontSize=\'24px\';" onmouseout="this.style.backgroundColor=\'transparent\'; this.style.fontSize=\'20px\';">{residue}</span>'
                            
                            # Add space every 10 residues for readability
                            if (i + 1) % 10 == 0:
                                seq_html += "  "
                            # Add line break every 100 residues
                            if (i + 1) % 100 == 0:
                                seq_html += "<br>"
                        seq_html += "</div>"
                        st.markdown(seq_html, unsafe_allow_html=True)
                        st.markdown("")
                else:
                    # Single chain protein
                    seq_html = "<div style='font-family: monospace; font-size: 20px; line-height: 1.2;'>"
                    for i, residue in enumerate(protein_seq):
                        # Calculate position (1-based)
                        pos = i + 1
                        # Add hover tooltip with residue code and index
                        seq_html += f'<span title="[ {residue} ] Index: {pos}" style="cursor: pointer; padding: 1px; border-radius: 2px; font-size: 20px;" onmouseover="this.style.backgroundColor=\'yellow\'; this.style.fontSize=\'24px\';" onmouseout="this.style.backgroundColor=\'transparent\'; this.style.fontSize=\'20px\';">{residue}</span>'
                        
                        # Add space every 10 residues for readability
                        if (i + 1) % 10 == 0:
                            seq_html += "  "
                        # Add line break every 50 residues
                        if (i + 1) % 100 == 0:
                            seq_html += "<br>"
                    seq_html += "</div>"
                    st.markdown(seq_html, unsafe_allow_html=True)
    else:
        st.info("No protein sequences available")

def parse_fasta_sequences(fasta_text: str) -> List[Tuple[str, str]]:
    """
    Parse FASTA format text and return list of (name, sequence) tuples.
    
    Args:
        fasta_text: FASTA format text with sequences
        
    Returns:
        List of tuples containing (sequence_name, sequence)
    """
    sequences = []
    current_name = None
    current_sequence = ""
    
    lines = fasta_text.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_name and current_sequence:
                sequences.append((current_name, current_sequence))
            
            # Start new sequence
            current_name = line[1:]  # Remove '>' prefix
            current_sequence = ""
        else:
            # Add to current sequence
            current_sequence += line
    
    # Add the last sequence
    if current_name and current_sequence:
        sequences.append((current_name, current_sequence))
    
    return sequences

def parse_smiles_list(smiles_text: str) -> List[Tuple[str, str]]:
    """
    Parse SMILES text in FASTA format and return list of (name, smiles) tuples.
    
    Args:
        smiles_text: FASTA format text with SMILES strings
        
    Returns:
        List of tuples containing (drug_name, smiles)
    """
    drugs = []
    current_name = None
    current_smiles = ""
    
    lines = smiles_text.strip().split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous drug if exists
            if current_name and current_smiles:
                drugs.append((current_name, current_smiles))
            
            # Start new drug
            current_name = line[1:]  # Remove '>' prefix
            current_smiles = ""
        else:
            # Add to current SMILES
            current_smiles += line
    
    # Add the last drug
    if current_name and current_smiles:
        drugs.append((current_name, current_smiles))
    
    return drugs

def validate_smiles(smiles: str) -> bool:
    """
    Basic SMILES validation using RDKit.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if valid, False otherwise
    """
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except ImportError:
        # If RDKit is not available, do basic string validation
        if not smiles or len(smiles) < 3:
            return False
        # Check for basic SMILES characters
        valid_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]{}@+-=#$%:./\\')
        return all(c in valid_chars for c in smiles)
    except Exception:
        return False

def validate_ccd_code(ccd_code: str) -> bool:
    """
    Basic validation for CCD codes.
    
    Args:
        ccd_code: CCD code string to validate
        
    Returns:
        True if valid format, False otherwise
    """
    if not ccd_code or not isinstance(ccd_code, str):
        return False
    
    # CCD codes are typically 3-letter codes (e.g., HEM, ATP, NAD)
    # They should be uppercase and contain only letters
    ccd_code = ccd_code.strip().upper()
    if len(ccd_code) < 2 or len(ccd_code) > 4:
        return False
    
    # Check if it contains only letters
    return ccd_code.isalpha()

def validate_protein_sequence(protein_seq: str) -> Tuple[bool, str, Dict[str, str], str]:
    """
    Validate protein sequence input for multi-chain format.
    Converts input to uppercase before validation and storage.
    Removes all whitespace (spaces, newlines, tabs) for robust parsing.
    Args:
        protein_seq (str): Protein sequence input
    Returns:
        tuple: (is_valid, error_message, chains_dict, upper_seq)
    """
    protein_seq = re.sub(r'\s+', '', protein_seq.upper())
    if not protein_seq.strip():
        return True, "", {}, protein_seq
    
    # Check for invalid characters (only uppercase letters and : allowed)
    invalid_chars = re.findall(r'[^A-Z:]', protein_seq)
    if invalid_chars:
        invalid_chars_str = ', '.join(set(invalid_chars))
        return False, f"Invalid characters found: {invalid_chars_str}. Only uppercase letters (A-Z) and colon (:) are allowed.", {}, protein_seq
    
    # Split by colon to get chains
    chains = protein_seq.split(':')
    
    # Check if we have too many chains (max 25 chains: A-W, Y-Z, X is reserved for ligand)
    if len(chains) > 25:
        return False, f"Too many protein chains ({len(chains)}). Maximum allowed is 25 chains (A-W, Y-Z). Chain X is reserved for ligands.", {}, protein_seq
    
    # Validate each chain contains only uppercase letters
    for i, chain in enumerate(chains):
        if not chain.strip():
            return False, f"Empty chain found at position {i+1}. Each chain must contain amino acid sequence.", {}, protein_seq
        if not re.match(r'^[A-Z]+$', chain.strip()):
            return False, f"Chain {i+1} contains invalid characters. Only uppercase letters are allowed.", {}, protein_seq
    
    # Create chains dictionary with incremental IDs (A-W, Y-Z, X reserved for ligand)
    chains_dict = {}
    for i, chain in enumerate(chains):
        if i < 23:  # A-W (0-22)
            chain_id = chr(65 + i)  # A, B, C, ..., W (65-87)
        else:  # Y-Z (23-24)
            chain_id = chr(89 + (i - 23))  # Y, Z (89-90)
        chains_dict[chain_id] = chain.strip()
    
    return True, "", chains_dict, protein_seq

def parse_manual_binding_pocket_residues(manual_residues: str) -> List[List[str]]:
    """
    Parse manual binding pocket residues input.
    
    Args:
        manual_residues: String like "A,123 B,45" with space-separated CHAIN_ID,RES_IDX pairs
        
    Returns:
        List of [chain_id, residue_idx] lists
    """
    contacts = []
    for entry in manual_residues.strip().split():
        if ',' in entry:
            parts = entry.split(',')
            if len(parts) == 2:
                contacts.append([parts[0].strip(), parts[1].strip()])
    return contacts

def find_sequence_matches_in_chains(sequence: str, chains_dict: Dict[str, str]) -> Tuple[List[str], Dict[str, List[int]]]:
    """
    Find which chains contain a given sequence and return matching positions.
    
    Args:
        sequence: Sequence to search for
        chains_dict: Dictionary mapping chain_id to sequence
        
    Returns:
        Tuple of (matching_chains, chain_to_match_indices)
    """
    matching_chains = []
    chain_to_match_indices = {}
    
    for chain_id, seq_str in chains_dict.items():
        if sequence:
            indices = []
            start = 0
            while True:
                idx = seq_str.find(sequence, start)
                if idx == -1:
                    break
                indices.append(idx + 1)  # 1-based
                start = idx + 1
            if indices:
                matching_chains.append(chain_id)
                chain_to_match_indices[chain_id] = indices
    
    return matching_chains, chain_to_match_indices

def generate_contacts_from_sequence_matches(sequence: str, selected_chains: List[str], 
                                          chain_to_match_indices: Dict[str, List[int]], 
                                          chains_dict: Dict[str, str]) -> List[List[str]]:
    """
    Generate contact list from sequence matches.
    
    Args:
        sequence: The sequence that was matched
        selected_chains: List of chain IDs to use
        chain_to_match_indices: Dictionary mapping chain_id to list of start positions
        chains_dict: Dictionary mapping chain_id to sequence
        
    Returns:
        List of [chain_id, residue_idx] contacts
    """
    contacts = []
    for chain in selected_chains:
        indices = chain_to_match_indices.get(chain, [])
        seq_str = chains_dict.get(chain, "")
        for idx in indices:
            if idx > 0 and idx + len(sequence) - 1 <= len(seq_str):
                for j in range(len(sequence)):
                    contacts.append([chain, str(idx + j)])
    return contacts

def display_binding_pocket_section(protein_sequence: str = None, default_constraints: Dict = None) -> Dict:
    """
    Display binding pocket constraints section and return the constraints dictionary.
    
    Args:
        protein_sequence: Current protein sequence (can be None or empty)
        default_constraints: Default constraints to use
        
    Returns:
        Dictionary with binding pocket constraints or None if not set
    """
    if default_constraints is None:
        default_constraints = {
            'binder': 'X',
            'contacts': [],
            'max_distance': 7.0,
            'mode': 'automated',
            'sequence': '',
            'manual_residues': ''
        }
    
    # Try to get the latest input from session state (Multi-Protein Mode)
    latest_protein_seq = None
    if 'protein_input_normal' in st.session_state and st.session_state['protein_input_normal']:
        latest_protein_seq = st.session_state['protein_input_normal']
    # Try to get the latest input from session state (Mutation Mode)
    elif 'wt_protein_input' in st.session_state and st.session_state['wt_protein_input']:
        latest_protein_seq = st.session_state['wt_protein_input']
    # Fallback to argument
    elif protein_sequence:
        latest_protein_seq = protein_sequence
    else:
        latest_protein_seq = ''

    st.markdown("""
    Specify the binding pocket residues for the ligand. This information will be used as constraints for structure prediction.
    """)

    col1, col2, _ = st.columns(3)
    with col1:
        pocket_mode = st.segmented_control(
            "Pocket constraint mode",
            options=["Automated (by sequence)", "Manual (by residue list)"],
            default="Automated (by sequence)",
            key="batch_pocket_mode"
        )

    with col2:
        binder_chain = "X"  # X is always the ligand
        max_distance = st.number_input(
            "Max distance (Å)",
            min_value=1.0,
            max_value=30.0,
            icon=":material/arrow_range:",
            value=float(default_constraints.get('max_distance', 7.0)),
            step=0.1,
            key="batch_max_distance_input"
    )

    # Parse protein sequence to get chains
    chains_dict = {}
    is_valid = False

    if latest_protein_seq and latest_protein_seq.strip():
        is_valid, _, chains_dict, _ = validate_protein_sequence(latest_protein_seq)

    contacts = []

    if pocket_mode == "Automated (by sequence)":
        st.markdown("**Automated Mode:** Enter one or more short sequences, search for them in your protein chains, and select the chain(s) to use for each.")

        # Initialize session state for sequence queries
        if 'batch_pocket_seq_queries' not in st.session_state:
            st.session_state['batch_pocket_seq_queries'] = [
                {'seq': default_constraints.get('sequence', ''), 'id': 0}
            ]

        seq_queries = st.session_state['batch_pocket_seq_queries']
        next_id = max([q['id'] for q in seq_queries], default=0) + 1
        all_contacts = []
        remove_ids = []

        for i, q in enumerate(seq_queries):
            cols = st.columns([3, 3, 1], vertical_alignment="bottom")
            with cols[0]:
                seq = st.text_input(f"Sequence to search for", value=q['seq'], key=f"batch_pocket_seq_query_{q['id']}")
                
                # Find matching chains for this sequence
                matching_chains, chain_to_match_indices = find_sequence_matches_in_chains(seq, chains_dict)
                
                with cols[1]:
                    selected_chains = st.multiselect(
                        "Select chain(s) containing the sequence",
                        options=matching_chains,
                        default=matching_chains,
                        key=f"batch_pocket_chain_multiselect_{q['id']}"
                    )
                
                with cols[2]:
                    trash_disabled = len(seq_queries) == 1
                    if st.button("", key=f"batch_pocket_seq_trash_{q['id']}", disabled=trash_disabled, icon=":material/delete:", type="tertiary"):
                        remove_ids.append(q['id'])
                
                # Prepare contacts for this sequence
                contacts = generate_contacts_from_sequence_matches(seq, selected_chains, chain_to_match_indices, chains_dict)
                
                if matching_chains:
                    all_residues = sorted(set([f'{c},{r}' for c, r in contacts]))
                    st.write(f"Residues matched: {all_residues}")
                elif seq:
                    st.info("No matches found for the sequence in any chain.")
                
                all_contacts.extend(contacts)
            
            # Remove any deleted queries
            if remove_ids:
                st.session_state['batch_pocket_seq_queries'] = [q for q in seq_queries if q['id'] not in remove_ids]
                st.rerun()
            
            # Add plus button
            if st.button("", key="batch_pocket_seq_add", icon=":material/add:", type="tertiary"):
                st.session_state['batch_pocket_seq_queries'].append({'seq': '', 'id': next_id})
                st.rerun()
            
            # Show all contacts summary
            if all_contacts:
                all_residues = sorted(set([f'{c},{r}' for c, r in all_contacts]))
                st.write(f"All residues matched: {all_residues}")
            
            _, col, _ = st.columns([2, 1, 2])
            with col:
                pocket_submitted = st.button("Apply Binding Pocket Constraints", key="batch_pocket_auto_submit", icon=":material/check_circle:")
            
            if pocket_submitted and all_contacts:
                st.success("Binding pocket constraints updated (automated mode).")
                return {
                    'binder': binder_chain.strip(),
                    'contacts': all_contacts,
                    'max_distance': float(max_distance),
                    'mode': 'automated',
                    'sequence': ','.join([q['seq'] for q in st.session_state['batch_pocket_seq_queries']])
                }
        
        else:  # Manual mode
            st.markdown("**Manual Mode:** Enter residues as `CHAIN_ID,RES_IDX` separated by spaces (e.g., `A,123 B,45`).")
            manual_residues = st.text_input(
                "Residues (space-separated)", 
                value=default_constraints.get('manual_residues', ''), 
                key="batch_pocket_manual_residues"
            )
            
            # Parse manual input into [[CHAIN_ID, RES_IDX], ...]
            contacts = parse_manual_binding_pocket_residues(manual_residues)
            
            _, col, _ = st.columns([2, 1, 2])
            with col:
                pocket_submitted = st.button("Apply Binding Pocket Constraints", key="batch_pocket_manual_submit", icon=":material/check_circle:")
            
            if pocket_submitted and contacts:
                st.success("Binding pocket constraints updated (manual mode).")
                return {
                    'binder': binder_chain.strip(),
                    'contacts': contacts,
                    'max_distance': float(max_distance),
                    'mode': 'manual',
                    'manual_residues': manual_residues
                }
    
    return None

def display_ptm_section(protein_sequence: str = None, default_ptms: Dict = None, protein_sequences: List[Tuple[str, str]] = None) -> Dict:
    """
    Display post-translational modifications section and return the PTM dictionary.
    
    Args:
        protein_sequence: Current protein sequence (can be None or empty)
        default_ptms: Default PTM modifications to use
        protein_sequences: List of (name, sequence) tuples from parsed input
        
    Returns:
        Dictionary with PTM modifications or None if not set
    """
    if default_ptms is None:
        default_ptms = {
            'modifications': []
        }
    
    # Try to get the latest input from session state (Multi-Protein Mode)
    latest_protein_seq = None
    if 'protein_input_normal' in st.session_state and st.session_state['protein_input_normal']:
        latest_protein_seq = st.session_state['protein_input_normal']
    # Try to get the latest input from session state (Mutation Mode)
    elif 'wt_protein_input' in st.session_state and st.session_state['wt_protein_input']:
        latest_protein_seq = st.session_state['wt_protein_input']
    # Fallback to argument
    elif protein_sequence:
        latest_protein_seq = protein_sequence
    else:
        latest_protein_seq = ''

    # Add sequence viewer button
    col_title, col_button = st.columns([3, 1], vertical_alignment="bottom")
    with col_title:
        st.markdown("Specify post-translational modifications (PTMs) for protein residues. These modifications will be included in the structure prediction. Click :material/visibility: View Sequences to see residue indices.")
    with col_button:
        if st.button("View Sequences", key="batch_ptm_view_sequences", icon=":material/visibility:", type="tertiary"):
            show_sequence_viewer_dialog(protein_sequences)
        
    # Initialize session state for PTM modifications
    if 'batch_ptm_modifications' not in st.session_state:
        st.session_state['batch_ptm_modifications'] = []
        
    ptm_modifications = st.session_state['batch_ptm_modifications']
    next_id = max([m['id'] for m in ptm_modifications], default=0) + 1
    remove_ids = []
        
    # Common PTM CCD codes with descriptions
    common_ptms = {
        "SEP": "Phosphoserine (phosphorylated serine)",
        "TPO": "Phosphothreonine (phosphorylated threonine)", 
        "PTR": "Phosphotyrosine (phosphorylated tyrosine)",
        "HID": "Phosphohistidine (phosphorylated histidine)",
        "ALY": "N6-acetyllysine (acetylated lysine)",
        "MLY": "N6-methyllysine (methylated lysine)",
        "MLZ": "N6,N6-dimethyllysine (dimethylated lysine)",
        "ML3": "N6,N6,N6-trimethyllysine (trimethylated lysine)",
        "ARG": "N-methylarginine (methylated arginine)",
        "MAL": "S-malonylcysteine (malonylated cysteine)",
        "HYP": "4-hydroxyproline (hydroxylated proline)",
        "HYK": "5-hydroxylysine (hydroxylated lysine)",
        "ASN": "N4-hydroxyasparagine (hydroxylated asparagine)",
        "PAL": "S-palmitoylcysteine (palmitoylated cysteine)",
        "SUC": "N-succinylasparagine (succinylated asparagine)",
        "SNO": "S-nitrosocysteine (S-nitrosylated cysteine)",
        "FOR": "N-formyltryptophan (formylated tryptophan)",
        "CRT": "N6-crotonyllysine (crotonylated lysine)",
        "CIT": "Citrulline (citrullinated arginine)"
    }
        
    # Input fields for new PTM
    col1, col2, col3, col4 = st.columns(4)
        
    with col1:
        # Protein selection dropdown
        if protein_sequences:
            protein_options = ["All models"] + [name for name, _ in protein_sequences]
            selected_protein = st.selectbox(
                "Apply to protein",
                options=protein_options,
                key="batch_ptm_protein_select",
                help="Select which protein model to apply the modification to"
            )
        else:
            selected_protein = "All models"
            st.write("**Apply to protein:** All models (no proteins loaded)")
        
    with col2:
        # Chain selection dropdown
        if selected_protein == "All models" and protein_sequences:
            # For "All models", collect all unique chains from all proteins
            all_chains = set()
            for protein_name, protein_seq in protein_sequences:
                is_valid_protein, _, chains_dict_protein, _ = validate_protein_sequence(protein_seq)
                if is_valid_protein:
                    if chains_dict_protein:
                        all_chains.update(chains_dict_protein.keys())
                    else:
                        all_chains.add('A')  # Single chain
                
            chain_options = ["All chains"] + sorted(list(all_chains))
            selected_chains = st.multiselect(
                "Apply to chains",
                options=chain_options,
                default=["All chains"],
                key="batch_ptm_chains_select",
                help="Select which chains to apply the modification to"
            )
        elif selected_protein != "All models" and protein_sequences:
            # For specific protein, get its chains
            chains_for_protein = []
            for protein_name, protein_seq in protein_sequences:
                if protein_name == selected_protein:
                    is_valid_protein, _, chains_dict_protein, _ = validate_protein_sequence(protein_seq)
                    if is_valid_protein:
                        if chains_dict_protein:
                            chains_for_protein = sorted(list(chains_dict_protein.keys()))
                        else:
                            chains_for_protein = ['A']  # Single chain
                    break
                
            chain_options = ["All chains"] + chains_for_protein
            selected_chains = st.multiselect(
                "Apply to chains",
                options=chain_options,
                default=["All chains"],
                key="batch_ptm_chains_select",
                help="Select which chains to apply the modification to"
            )
        else:
            selected_chains = ["All chains"]
            st.write("**Apply to chains:** All chains")
        
    with col3:
        # Index input
        position = st.number_input(
            "Residue index (starting at 1)",
            min_value=1,
            value=1,
            key="batch_ptm_position_input",
            help="Enter the residue index (1-based indexing). Click View Sequences in the upper right corner to identify the index. This is different from residue number."
        )
        
    with col4:
        # CCD code selection (dropdown only, no search)
        ccd_options = list(common_ptms.keys())
        ccd_options.sort()
            
        selected_ccd = st.selectbox(
            "PTM type (CCD code)",
            options=ccd_options,
            format_func=lambda x: f"{x} - {common_ptms.get(x, 'Custom')}",
            key="batch_ptm_ccd_select"
        )

    # Validation and context display
    verification_data = []  # Initialize verification_data
        
    # Check if we have protein sequences to validate against
    if protein_sequences and len(protein_sequences) > 0:
        # Create verification table for the specified position
            
        # Determine which chains to check based on selection
        chains_to_check = []
        if "All chains" in selected_chains:
            chains_to_check = ["All chains"]
        else:
            chains_to_check = selected_chains
            
        if selected_protein == "All models":
            # Check all protein models
            for protein_name, protein_seq in protein_sequences:
                is_valid_protein, _, chains_dict_protein, _ = validate_protein_sequence(protein_seq)
                if is_valid_protein:
                    if chains_dict_protein:
                        # Multi-chain protein
                        for chain_id in sorted(chains_dict_protein.keys()):
                            # Check if this chain should be included
                            if "All chains" in chains_to_check or chain_id in chains_to_check:
                                chain_seq = chains_dict_protein[chain_id]
                                if 1 <= position <= len(chain_seq):
                                    residue = chain_seq[position - 1]
                                    context = create_residue_context(chain_seq, position - 1, residue)
                                    verification_data.append({
                                        "Protein": protein_name,
                                        "Chain": chain_id,
                                        "Index": position,
                                        "Residue": residue,
                                        "Context": context,
                                        "Status": "Valid"
                                    })
                                else:
                                    verification_data.append({
                                        "Protein": protein_name,
                                        "Chain": chain_id,
                                        "Index": position,
                                        "Residue": "Out of range",
                                        "Context": f"Index {position} out of range (1-{len(chain_seq)})",
                                        "Status": "Invalid"
                                    })
                    else:
                        # Single chain protein
                        if "All chains" in chains_to_check or 'A' in chains_to_check:
                            if 1 <= position <= len(protein_seq):
                                residue = protein_seq[position - 1]
                                context = create_residue_context(protein_seq, position - 1, residue)
                                verification_data.append({
                                    "Protein": protein_name,
                                    "Chain": "A",
                                    "Index": position,
                                    "Residue": residue,
                                    "Context": context,
                                    "Status": "Valid"
                                })
                            else:
                                verification_data.append({
                                    "Protein": protein_name,
                                    "Chain": "A",
                                    "Index": position,
                                    "Residue": "Out of range",
                                    "Context": f"Index {position} out of range (1-{len(protein_seq)})",
                                    "Status": "Invalid"
                                })
        else:
            # Check specific protein model
            for protein_name, protein_seq in protein_sequences:
                if protein_name == selected_protein:
                    is_valid_protein, _, chains_dict_protein, _ = validate_protein_sequence(protein_seq)
                    if is_valid_protein:
                        if chains_dict_protein:
                            # Multi-chain protein
                            for chain_id in sorted(chains_dict_protein.keys()):
                                # Check if this chain should be included
                                if "All chains" in chains_to_check or chain_id in chains_to_check:
                                    chain_seq = chains_dict_protein[chain_id]
                                    if 1 <= position <= len(chain_seq):
                                        residue = chain_seq[position - 1]
                                        context = create_residue_context(chain_seq, position - 1, residue)
                                        verification_data.append({
                                            "Protein": protein_name,
                                            "Chain": chain_id,
                                            "Index": position,
                                            "Residue": residue,
                                            "Context": context,
                                            "Status": "Valid"
                                        })
                                    else:
                                        verification_data.append({
                                            "Protein": protein_name,
                                            "Chain": chain_id,
                                            "Index": position,
                                            "Residue": "Out of range",
                                            "Context": f"Index {position} out of range (1-{len(chain_seq)})",
                                            "Status": "Invalid"
                                        })
                        else:
                            # Single chain protein
                            if "All chains" in chains_to_check or 'A' in chains_to_check:
                                if 1 <= position <= len(protein_seq):
                                    residue = protein_seq[position - 1]
                                    context = create_residue_context(protein_seq, position - 1, residue)
                                    verification_data.append({
                                        "Protein": protein_name,
                                        "Chain": "A",
                                        "Index": position,
                                        "Residue": residue,
                                        "Context": context,
                                        "Status": "Valid"
                                    })
                                else:
                                    verification_data.append({
                                        "Protein": protein_name,
                                        "Chain": "A",
                                        "Index": position,
                                        "Residue": "Out of range",
                                        "Context": f"Index {position} out of range (1-{len(protein_seq)})",
                                        "Status": "Invalid"
                                    })
                    break
            
        # Display verification table and current PTM modifications side by side
        if verification_data:        
            # Display tables side by side
            col1, col2 = st.columns([1, 1])
            with col1:
                st.markdown("##### :material/search: Selected Residue Indices")
                
                # Create verification DataFrame
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
                        "Protein": st.column_config.TextColumn("Protein", width="medium"),
                        "Chain": st.column_config.TextColumn("Chain", width="small"),
                        "Index": st.column_config.NumberColumn("Index", width="small"),
                        "Residue": st.column_config.TextColumn("Residue", width="small"),
                        "Context": st.column_config.TextColumn("Context", width="medium", help="Shows nearby residues for context"),
                        "Status": st.column_config.TextColumn("Status", width="small")
                    }
                )
                
            with col2:
                st.markdown("##### :material/verified: Current PTM Modifications")
                if ptm_modifications:
                    # Display current PTM modifications as simple text with inline delete buttons
                    for mod in ptm_modifications:
                        chain_id = mod.get('chain_id', 'A')
                        position = mod.get('position', '')
                        ccd = mod.get('ccd', '')
                        protein = mod.get('protein', 'All models')
                        description = common_ptms.get(ccd, f"Custom modification: {ccd}")
                            
                        # Create inline display with delete button
                        col_text, col_delete = st.columns([4, 1])
                        with col_text:
                            st.write(f"**{protein}**, Chain {chain_id}, Index {position}: {ccd} - {description}")
                        with col_delete:
                            if st.button("", key=f"batch_ptm_delete_{mod['id']}", icon=":material/delete:", type="tertiary"):
                                remove_ids.append(mod['id'])
                else:
                    st.info("No PTM modifications added yet")
                
            # Remove deleted modifications
            if remove_ids:
                st.session_state['batch_ptm_modifications'] = [m for m in ptm_modifications if m['id'] not in remove_ids]
                st.rerun()
        
    # Show summary
    valid_count = len([r for r in verification_data if r['Status'] == 'Valid'])
    invalid_count = len([r for r in verification_data if r['Status'] == 'Invalid'])
        
    if invalid_count > 0:
        st.warning(f"{invalid_count} index(es) are out of range. Please check your input.")
    else:
        if valid_count > 0:
            st.success(f"All {valid_count} index(es) are valid! Click 'Add PTM modification' to add to list of modifications to consider. Once done, click 'Apply PTM Modifications' to apply changes.")
        
    # Add button
    col1, col2, col3 = st.columns([1, 1, 1])
    with col2:
        if st.button("Add PTM modification", key="batch_ptm_add", icon=":material/add:", type="tertiary"):
            # Validate input
            if selected_ccd and position > 0 and selected_chains:
                # Handle multiple chains
                chains_to_add = []
                if "All chains" in selected_chains:
                    # For "All chains", we need to determine which chains actually exist
                    if selected_protein == "All models":
                        # Collect all chains from all proteins
                        all_chains = set()
                        for protein_name, protein_seq in protein_sequences:
                            is_valid_protein, _, chains_dict_protein, _ = validate_protein_sequence(protein_seq)
                            if is_valid_protein:
                                if chains_dict_protein:
                                    all_chains.update(chains_dict_protein.keys())
                                else:
                                    all_chains.add('A')
                        chains_to_add = sorted(list(all_chains))
                    else:
                        # Get chains for specific protein
                        for protein_name, protein_seq in protein_sequences:
                            if protein_name == selected_protein:
                                is_valid_protein, _, chains_dict_protein, _ = validate_protein_sequence(protein_seq)
                                if is_valid_protein:
                                    if chains_dict_protein:
                                        chains_to_add = sorted(list(chains_dict_protein.keys()))
                                    else:
                                        chains_to_add = ['A']
                                break
                else:
                    chains_to_add = selected_chains
                    
                # Check for existing modifications
                existing_modifications = []
                for chain_id in chains_to_add:
                    existing = [m for m in ptm_modifications 
                              if m.get('protein') == selected_protein and
                              m.get('chain_id') == chain_id and 
                              m.get('position') == position]
                    existing_modifications.extend(existing)
                    
                if existing_modifications:
                    existing_chains = [m.get('chain_id') for m in existing_modifications]
                    st.error(f"PTM modification at index {position} for {selected_protein}, chains {existing_chains} already exists!")
                else:
                    # Add new modifications for each chain
                    added_count = 0
                    for chain_id in chains_to_add:
                        new_mod = {
                            'id': next_id + added_count,
                            'chain_id': chain_id,
                            'position': position,
                            'ccd': selected_ccd,
                            'protein': selected_protein
                        }
                        st.session_state['batch_ptm_modifications'].append(new_mod)
                        added_count += 1
                        
                    chain_text = "all chains" if len(chains_to_add) > 1 else f"chain {chains_to_add[0]}"
                    st.success(f"Added {selected_ccd} modification at index {position} for {selected_protein}, {chain_text}")
                    st.rerun()
            else:
                st.error("Please select a valid PTM type, index, and chains")

    # Apply button
    if ptm_modifications:
        _, col, _ = st.columns([2, 1, 2])
        with col:
            ptm_submitted = st.button("Apply PTM Modifications", key="batch_ptm_submit", icon=":material/check_circle:")
            
        if ptm_submitted:
            st.success(f"PTM modifications applied: {len(ptm_modifications)} modification(s)")
            return {
                'modifications': ptm_modifications
            }

    return None

def display_protein_drug_filter_section(protein_sequences: List[Tuple[str, str]] = None, drug_sequences: List[Tuple[str, str]] = None) -> Dict:
    """
    Display protein-drug pairing filter section and return the filter dictionary.
    
    Args:
        protein_sequences: List of (name, sequence) tuples from parsed protein input
        drug_sequences: List of (name, smiles) tuples from parsed drug input
        
    Returns:
        Dictionary with protein-drug pair filters or None if not set
    """
    
    st.markdown("Specify which protein-drug combinations should be evaluated. If this field is filled out, only the specified protein-drug pairs will be evaluated, skipping all other combinations.")
        
    # Initialize session state for protein-drug filter
    if 'batch_protein_drug_filter' not in st.session_state:
        st.session_state['batch_protein_drug_filter'] = ""
        
    # Text input for protein-drug pairs
    filter_text = st.text_area(
        "Protein-drug pairs",
        value=st.session_state['batch_protein_drug_filter'],
        height=100,
        placeholder="WT, Dofetilide; G648S, Moxifloxacin; Y652A, E4031",
        help="Format: protein_name, drug_name; protein_name2, drug_name2\nYou can use semicolons (;) or newlines to separate pairs.\nSpaces around commas and semicolons will be stripped.",
        key="batch_protein_drug_filter_input"
    )
        
    # Update session state
    st.session_state['batch_protein_drug_filter'] = filter_text
        
    # Parse and validate the input
    parsed_pairs = []
    validation_errors = []
        
    if filter_text.strip():
        # Split by semicolons or newlines
        pairs_text = filter_text.strip()
        pairs_text = pairs_text.replace('\n', ';')  # Convert newlines to semicolons
        pair_strings = [pair.strip() for pair in pairs_text.split(';') if pair.strip()]
            
        available_proteins = set()
        available_drugs = set()
            
        # Get available protein and drug names
        if protein_sequences:
            available_proteins = {name for name, _ in protein_sequences}
        if drug_sequences:
            available_drugs = {name for name, _ in drug_sequences}
            
        for pair_str in pair_strings:
            if ',' in pair_str:
                parts = [part.strip() for part in pair_str.split(',')]
                if len(parts) == 2:
                    protein_name, drug_name = parts
                    parsed_pairs.append((protein_name, drug_name))
                        
                    # Validate protein and drug names
                    if protein_sequences and protein_name not in available_proteins:
                        validation_errors.append(f"Protein '{protein_name}' not found in protein sequences")
                    if drug_sequences and drug_name not in available_drugs:
                        validation_errors.append(f"Drug '{drug_name}' not found in drug sequences")
                else:
                    validation_errors.append(f"Invalid format: '{pair_str}' - expected 'protein, drug'")
            else:
                validation_errors.append(f"Invalid format: '{pair_str}' - missing comma separator")
        
    # Display validation results
    if filter_text.strip():
        if validation_errors:
            st.markdown("##### :material/error: Validation Errors")
            for error in validation_errors:
                st.error(error)
        elif parsed_pairs:
            st.success("All pairs are valid!")
        else:
            st.info("No valid pairs found")
        
    # Return the filter data if there are valid pairs
    if parsed_pairs and not validation_errors:
        return {
            'enabled': True,
            'pairs': parsed_pairs
        }
    elif filter_text.strip() and not parsed_pairs:
        return {
            'enabled': True,
            'pairs': []
        }

    return None

def calculate_filtered_job_count(protein_sequences: List[Tuple[str, str]], drug_sequences: List[Tuple[str, str]] = None, 
                                structure_only: bool = False, protein_drug_filter: Dict = None) -> int:
    """
    Calculate the total number of jobs that will be executed considering the protein-drug filter.
    
    Args:
        protein_sequences: List of (name, sequence) tuples
        drug_sequences: List of (name, smiles) tuples (None for structure-only mode)
        structure_only: Whether running in structure-only mode
        protein_drug_filter: Filter dictionary from display_protein_drug_filter_section
        
    Returns:
        Number of jobs that will be executed
    """
    if structure_only:
        # In structure-only mode, count proteins that pass the filter
        if not protein_drug_filter or not protein_drug_filter.get('enabled'):
            return len(protein_sequences)
        
        filter_pairs = protein_drug_filter.get('pairs', [])
        if not filter_pairs:
            return 0
        
        # Count proteins that appear in any filter pair
        filtered_proteins = set(pair[0] for pair in filter_pairs)
        return sum(1 for protein_name, _ in protein_sequences if protein_name in filtered_proteins)
    else:
        # In multi-protein mode, count protein-drug combinations that pass the filter
        if not protein_drug_filter or not protein_drug_filter.get('enabled'):
            return len(protein_sequences) * len(drug_sequences)
        
        filter_pairs = protein_drug_filter.get('pairs', [])
        if not filter_pairs:
            return 0
        
        # Count valid protein-drug combinations
        filter_set = set(filter_pairs)
        count = 0
        for protein_name, _ in protein_sequences:
            for drug_name, _ in drug_sequences:
                if (protein_name, drug_name) in filter_set:
                    count += 1
        return count

def should_evaluate_protein_drug_pair(protein_name: str, drug_name: str, protein_drug_filter: Dict = None) -> bool:
    """
    Check if a protein-drug pair should be evaluated based on the filter.
    
    Args:
        protein_name: Name of the protein
        drug_name: Name of the drug (can be None for structure-only mode)
        protein_drug_filter: Filter dictionary from display_protein_drug_filter_section
        
    Returns:
        True if the pair should be evaluated, False if it should be skipped
    """
    # If no filter is set, evaluate all pairs
    if not protein_drug_filter or not protein_drug_filter.get('enabled'):
        return True
    
    # If filter is enabled but no pairs are specified, skip all
    filter_pairs = protein_drug_filter.get('pairs', [])
    if not filter_pairs:
        return False
    
    # For structure-only mode (no drug), check if protein is in any pair
    if drug_name is None:
        return any(pair[0] == protein_name for pair in filter_pairs)
    
    # For multi-protein mode, check if the exact protein-drug pair exists
    return (protein_name, drug_name) in filter_pairs

def parse_mutations(mutation_text: str) -> List[List[Tuple[str, int, str]]]:
    """
    Parse mutation text and return list of mutation lists.
    
    Args:
        mutation_text: Mutation text like "A102S,G99E" where first letter is wild-type residue
        
    Returns:
        List of mutation lists, where each inner list represents one mutant
        Each mutation is a tuple of (chain_id, residue_number, new_residue)
    """
    if not mutation_text.strip():
        return []
    
    mutants = []
    # Split by comma to get individual mutants
    mutant_strings = [s.strip() for s in mutation_text.split(',')]
    
    for mutant_str in mutant_strings:
        if not mutant_str:
            continue
            
        mutations = []
        # Split by hyphen to get multiple mutations in one mutant
        mutation_parts = [s.strip() for s in mutant_str.split('-')]
        
        for part in mutation_parts:
            if not part or len(part) < 3:
                continue
                
            # Parse format like "A102S" where A is wild-type, 102 is position, S is new residue
            # Extract wild-type residue (first character)
            wt_residue = part[0].upper()
            # Extract residue number (middle part)
            residue_part = part[1:-1]
            # Extract new residue (last character)
            new_residue = part[-1].upper()
            
            # Validate wild-type residue (should be a valid amino acid)
            valid_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
            if wt_residue not in valid_amino_acids or new_residue not in valid_amino_acids:
                continue
            
            try:
                residue_number = int(residue_part)
                if residue_number > 0:  # Ensure positive residue number
                    # For now, assume single chain (chain A) - could be extended for multi-chain
                    chain_id = 'A'
                    mutations.append((chain_id, residue_number, new_residue))
            except ValueError:
                # Skip invalid mutations
                continue
        
        if mutations:
            mutants.append(mutations)
    
    return mutants

def apply_mutations_to_sequence(wt_sequence: str, mutations: List[Tuple[str, int, str]], chain_starts: Dict[str, int], chains_dict: Dict[str, str] = None) -> str:
    """
    Apply mutations to wild-type sequence.
    
    Args:
        wt_sequence: Wild-type protein sequence
        mutations: List of (chain_id, residue_number, new_residue) tuples
        chain_starts: Dictionary mapping chain_id to starting residue number
        chains_dict: Dictionary mapping chain_id to sequence (for multi-chain proteins)
        
    Returns:
        Mutated sequence (maintains multi-chain format with colons if original was multi-chain)
    """
    # For multi-chain proteins, work with individual chains
    if chains_dict and len(chains_dict) > 1:
        # Create a copy of chains_dict to modify
        mutated_chains = {chain_id: list(chain_seq) for chain_id, chain_seq in chains_dict.items()}
        
        for chain_id, residue_number, new_residue in mutations:
            # For multi-chain proteins, we need to handle both heteromultimers and homomultimers
            # In homomultimers (e.g., homotetramers), all chains start at the same residue number
            # and mutations should be applied to all chains that have the specified residue number
            
            # Find all chains that contain this residue number
            target_chains = []
            if chain_id == 'A':  # Default chain, need to detect which chains to target
                for cid in sorted(chains_dict.keys()):
                    if cid in chain_starts:
                        chain_start = chain_starts[cid]
                        chain_length = len(chains_dict[cid])
                        chain_end = chain_start + chain_length - 1
                        
                        if chain_start <= residue_number <= chain_end:
                            target_chains.append(cid)
            else:
                # Specific chain was specified, only target that chain
                target_chains = [chain_id]
            
            # Apply mutation to all target chains
            for target_chain_id in target_chains:
                if target_chain_id in chain_starts and target_chain_id in mutated_chains:
                    # Calculate position within the specific chain
                    chain_start = chain_starts[target_chain_id]
                    position_in_chain = residue_number - chain_start
                    
                    # Check if position is valid within the chain
                    if 0 <= position_in_chain < len(mutated_chains[target_chain_id]):
                        # Validate that the new residue is a valid amino acid
                        if new_residue.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                            mutated_chains[target_chain_id][position_in_chain] = new_residue.upper()
        
        # Reconstruct the multi-chain sequence with colons
        # We need to join the chains in the same order as the original input
        # Since we don't have the original order, we'll use the order in chains_dict
        result_parts = []
        for chain_id in chains_dict.keys():
            if chain_id in mutated_chains:
                result_parts.append(''.join(mutated_chains[chain_id]))
        
        # Join with colons to maintain multi-chain format
        return ':'.join(result_parts)
    else:
        # Single chain or fallback - use the original logic
        seq_list = list(wt_sequence)
        
        for chain_id, residue_number, new_residue in mutations:
            if chain_id in chain_starts:
                # Calculate position in concatenated sequence
                chain_start = chain_starts[chain_id]
                absolute_pos = residue_number - chain_start
                
                # Check if position is valid
                if 0 <= absolute_pos < len(seq_list):
                    # Validate that the new residue is a valid amino acid
                    if new_residue.upper() in 'ACDEFGHIKLMNPQRSTVWY':
                        seq_list[absolute_pos] = new_residue.upper()
        
        return ''.join(seq_list)

def generate_mutant_name_from_text(mutation_text: str) -> str:
    """
    Generate a name for the mutant from the original mutation text.
    
    Args:
        mutation_text: Original mutation text like "A102S,G99E"
        
    Returns:
        Mutant name string
    """
    if not mutation_text.strip():
        return "WT"
    
    # Parse the mutation text to extract the format we want for the name
    mutation_strings = []
    mutant_strings = [s.strip() for s in mutation_text.split(',')]
    
    for mutant_str in mutant_strings:
        if not mutant_str:
            continue
            
        # Split by - to get multiple mutations in one mutant
        mutation_parts = [s.strip() for s in mutant_str.split('-')]
        
        for part in mutation_parts:
            if not part or len(part) < 3:
                continue
                
            # Extract wild-type residue, residue number, and new residue
            wt_residue = part[0].upper()
            residue_part = part[1:-1]
            new_residue = part[-1].upper()
            
            try:
                residue_number = int(residue_part)
                if residue_number > 0:
                    mutation_strings.append(f"{wt_residue}{residue_number}{new_residue}")
            except ValueError:
                continue
    
    return "-".join(mutation_strings) if mutation_strings else "WT"

def verify_mutations_with_wt_residues(wt_sequence: str, mutation_text: str, chain_starts: Dict[str, int], chains_dict: Dict[str, str] = None) -> List[Tuple[bool, str, str, str, int, str, str]]:
    """
    Verify mutations by parsing mutation text that includes wild-type residues.
    Expected format: "A102S,G99E" where A and G are the wild-type residues.
    
    Args:
        wt_sequence: Wild-type protein sequence
        mutation_text: Mutation text like "A102S,G99E" where first letter is wild-type residue
        chain_starts: Dictionary mapping chain_id to starting residue number
        chains_dict: Dictionary mapping chain_id to sequence (for multi-chain proteins)
        
    Returns:
        List of (is_valid, chain_id, residue_number, wt_residue, actual_residue, new_residue, context) tuples
    """
    verification_results = []
    
    if not mutation_text.strip():
        return verification_results
    
    # Split by comma to get individual mutations
    mutation_strings = [s.strip() for s in mutation_text.split(',')]
    
    for mutation_str in mutation_strings:
        if not mutation_str:
            continue
            
        # Split by hyphen to get multiple mutations in one mutant
        mutation_parts = [s.strip() for s in mutation_str.split('-')]
        
        for part in mutation_parts:
            if not part or len(part) < 3: 
                continue
                
            # Parse format like "A102S" where A is wild-type, 102 is position, S is new residue
            # Extract wild-type residue (first character)
            wt_residue = part[0].upper()
            # Extract residue number (middle part)
            residue_part = part[1:-1]
            # Extract new residue (last character)
            new_residue = part[-1].upper()
            
            # Validate wild-type residue (should be a valid amino acid)
            valid_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
            if wt_residue not in valid_amino_acids or new_residue not in valid_amino_acids:
                continue
            
            try:
                residue_number = int(residue_part)
                if residue_number <= 0:
                    continue
                    
                # For multi-chain proteins, check all chains that contain this residue number
                # This handles both heteromultimers and homomultimers
                if chains_dict and len(chains_dict) > 1:
                    # Find all chains that contain this residue number
                    target_chains = []
                    for chain_id in sorted(chains_dict.keys()):
                        if chain_id in chain_starts:
                            chain_start = chain_starts[chain_id]
                            chain_length = len(chains_dict[chain_id])
                            chain_end = chain_start + chain_length - 1
                            
                            if chain_start <= residue_number <= chain_end:
                                target_chains.append(chain_id)
                    
                    # Verify mutation against all target chains
                    for chain_id in target_chains:
                        chain_start = chain_starts[chain_id]
                        chain_sequence = chains_dict[chain_id]
                        position_in_chain = residue_number - chain_start
                        
                        # Check if position is valid within the chain
                        if 0 <= position_in_chain < len(chain_sequence):
                            actual_residue = chain_sequence[position_in_chain]
                            
                            # Create context string with nearby residues from the chain
                            context = create_residue_context(chain_sequence, position_in_chain, actual_residue)
                            
                            # Check if the wild-type residue matches what's in the sequence
                            is_valid = (actual_residue == wt_residue)
                            verification_results.append((is_valid, chain_id, residue_number, wt_residue, actual_residue, new_residue, context))
                        else:
                            # Position out of range for this chain
                            context = f"Position {residue_number} out of range for chain {chain_id} (range: {chain_start}-{chain_start + len(chain_sequence) - 1})"
                            verification_results.append((False, chain_id, residue_number, wt_residue, "Position out of range", new_residue, context))
                else:
                    # Single chain protein - use original logic
                    chain_id = 'A'  # Default
                    if chain_id in chain_starts and chain_id in chains_dict:
                        # Calculate position within the specific chain
                        chain_start = chain_starts[chain_id]
                        chain_sequence = chains_dict[chain_id]
                        position_in_chain = residue_number - chain_start
                        
                        # Check if position is valid within the chain
                        if 0 <= position_in_chain < len(chain_sequence):
                            actual_residue = chain_sequence[position_in_chain]
                            
                            # Create context string with nearby residues from the chain
                            context = create_residue_context(chain_sequence, position_in_chain, actual_residue)
                            
                            # Check if the wild-type residue matches what's in the sequence
                            is_valid = (actual_residue == wt_residue)
                            verification_results.append((is_valid, chain_id, residue_number, wt_residue, actual_residue, new_residue, context))
                        else:
                            # Position out of range for this chain
                            context = f"Position {residue_number} out of range for chain {chain_id} (range: {chain_start}-{chain_start + len(chain_sequence) - 1})"
                            verification_results.append((False, chain_id, residue_number, wt_residue, "Position out of range", new_residue, context))
                    else:
                        # Chain not found
                        context = f"Chain {chain_id} not found"
                        verification_results.append((False, chain_id, residue_number, wt_residue, f"Chain {chain_id} not found", new_residue, context))
                    
            except ValueError:
                # Skip invalid mutations
                continue
    
    return verification_results

def create_residue_context(sequence: str, position: int, residue: str, context_size: int = 2) -> str:
    """
    Create a context string showing nearby residues around the specified position.
    
    Args:
        sequence: Protein sequence
        position: Position of the residue (0-based)
        residue: The residue at the position
        context_size: Number of residues to show on each side
        
    Returns:
        Context string like "ACK→S←RFV" where S is the residue at the position
    """
    if position < 0 or position >= len(sequence):
        return residue
    
    # Calculate start and end positions for context
    start_pos = max(0, position - context_size)
    end_pos = min(len(sequence), position + context_size + 1)
    
    # Get the context sequence
    context_seq = sequence[start_pos:end_pos]
    
    # Calculate the position of the residue within the context
    residue_pos_in_context = position - start_pos
    
    # Build the context string
    if residue_pos_in_context == 0:
        # Residue is at the beginning of the context
        context_str = f"{residue} ← {context_seq[1:]}"
    elif residue_pos_in_context == len(context_seq) - 1:
        # Residue is at the end of the context
        context_str = f"{context_seq[:-1]} → {residue}"
    else:
        # Residue is in the middle
        before = context_seq[:residue_pos_in_context]
        after = context_seq[residue_pos_in_context + 1:]
        context_str = f"{before} → {residue} ← {after}"
    
    return context_str


def calculate_position_in_chain(residue_number: int, chain_id: str, chain_starts: Dict[str, int]) -> int:
    """
    Calculate the 0-based position within a specific chain for a given residue number.
    
    Args:
        residue_number: The residue number
        chain_id: The chain ID
        chain_starts: Dictionary mapping chain_id to starting residue number
        
    Returns:
        0-based position within the chain
    """
    if chain_id in chain_starts:
        chain_start = chain_starts[chain_id]
        return residue_number - chain_start
    return -1


# Mutation Querying and Display Functions

def get_protein_info(sequence: str) -> Dict:
    """
    Identify protein using multiple search strategies via UniProt API.
    Returns basic protein information for mutation lookup.
    """
    try:
        # Clean sequence
        clean_seq = re.sub(r'[^A-Z]', '', sequence.upper())
        if len(clean_seq) < 20:
            return {"error": "Sequence too short for reliable identification"}

        seq_length = len(clean_seq)

        # Strategy 1: Direct sequence search is no longer supported by UniProt API
        # Skip to keyword and fragment-based searches
        url = "https://rest.uniprot.org/uniprotkb/search"

        # Strategy 2: Search by protein family keywords for common proteins
        # First, try to identify protein type from sequence patterns
        kinase_patterns = ['GSGAFG', 'VAIK', 'DFG', 'APE']  # Common kinase motifs
        is_likely_kinase = any(motif in clean_seq for motif in kinase_patterns)
        
        if is_likely_kinase:
            common_protein_patterns = [
                ('kinase', ['kinase', 'protein kinase', 'tyrosine kinase', 'serine kinase', 'MAP kinase']),
                ('hemoglobin', ['hemoglobin', 'alpha', 'beta']),
                ('myoglobin', ['myoglobin']),
                ('insulin', ['insulin']),
                ('lysozyme', ['lysozyme']),
                ('trypsin', ['trypsin']),
            ]
        else:
            common_protein_patterns = [
                ('hemoglobin', ['hemoglobin', 'alpha', 'beta']),
                ('myoglobin', ['myoglobin']),
                ('insulin', ['insulin']),
                ('lysozyme', ['lysozyme']),
                ('trypsin', ['trypsin']),
                ('kinase', ['kinase', 'protein kinase']),
            ]

        for protein_type, keywords in common_protein_patterns:
            for keyword in keywords:
                keyword_params = {
                    'query': f'reviewed:true AND {keyword} AND length:[{seq_length-20} TO {seq_length+20}]',
                    'format': 'json',
                    'size': 50,  # Increased to check more candidates
                    'fields': 'accession,id,protein_name,organism_name,length,sequence'
                }

                response = requests.get(url, params=keyword_params, timeout=15)
                if response.status_code == 200:
                    data = response.json()
                    results = data.get('results', [])

                    best_match = None
                    best_similarity = 0

                    for result in results:
                        seq = result.get('sequence', {}).get('value', '')
                        if seq and len(seq) > 0:
                            # Use a better similarity calculation
                            shorter_len = min(len(clean_seq), len(seq))
                            longer_len = max(len(clean_seq), len(seq))

                            # Calculate alignment-like similarity
                            matches = 0
                            for i in range(shorter_len):
                                if i < len(clean_seq) and i < len(seq) and clean_seq[i] == seq[i]:
                                    matches += 1

                            # Penalize length differences
                            length_penalty = abs(len(clean_seq) - len(seq)) / longer_len
                            similarity = (matches / shorter_len) * (1 - length_penalty * 0.5)

                            if similarity > best_similarity and similarity > 0.8:  # Higher threshold for better matches
                                best_similarity = similarity
                                best_match = result

                    if best_match:
                        return {
                            "accession": best_match.get('primaryAccession', ''),
                            "protein_name": best_match.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown'),
                            "organism": best_match.get('organism', {}).get('scientificName', 'Unknown'),
                            "similarity": best_similarity,
                            "length": best_match.get('sequence', {}).get('length', 0)
                        }

        # Strategy 3: Broad search with sequence fragments
        # Use first and last 20 amino acids as signature
        if len(clean_seq) >= 40:
            n_terminus = clean_seq[:20]
            c_terminus = clean_seq[-20:]

            fragment_params = {
                'query': f'reviewed:true AND length:[{seq_length-50} TO {seq_length+50}]',
                'format': 'json',
                'size': 100,
                'fields': 'accession,id,protein_name,organism_name,length,sequence'
            }

            response = requests.get(url, params=fragment_params, timeout=20)
            if response.status_code == 200:
                data = response.json()
                results = data.get('results', [])

                best_match = None
                best_similarity = 0

                for result in results:
                    seq = result.get('sequence', {}).get('value', '')
                    if seq and len(seq) >= 40:
                        # Check N and C terminal similarity
                        result_n_terminus = seq[:20]
                        result_c_terminus = seq[-20:]

                        n_matches = sum(1 for i in range(20) if n_terminus[i] == result_n_terminus[i])
                        c_matches = sum(1 for i in range(20) if c_terminus[i] == result_c_terminus[i])

                        terminal_similarity = (n_matches + c_matches) / 40.0

                        # Also check internal regions
                        internal_matches = 0
                        internal_total = 0
                        step = max(1, len(clean_seq) // 20)  # Sample every step residues

                        for i in range(20, min(len(clean_seq), len(seq)) - 20, step):
                            if clean_seq[i] == seq[i]:
                                internal_matches += 1
                            internal_total += 1

                        internal_similarity = internal_matches / max(1, internal_total)

                        # Combined similarity score
                        combined_similarity = terminal_similarity * 0.6 + internal_similarity * 0.4

                        if combined_similarity > best_similarity and combined_similarity > 0.3:
                            best_similarity = combined_similarity
                            best_match = result

                if best_match:
                    return {
                        "accession": best_match.get('primaryAccession', ''),
                        "protein_name": best_match.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown'),
                        "organism": best_match.get('organism', {}).get('scientificName', 'Unknown'),
                        "similarity": best_similarity,
                        "length": best_match.get('sequence', {}).get('length', 0)
                    }

        # Strategy 4: If still no match, try a very broad search and check each result
        if is_likely_kinase:
            broad_params = {
                'query': f'reviewed:true AND kinase AND length:[{seq_length-30} TO {seq_length+30}]',
                'format': 'json',
                'size': 200,  # Check many results
                'fields': 'accession,protein_name,organism_name,length,sequence'
            }
            
            try:
                response = requests.get(url, params=broad_params, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    results = data.get('results', [])
                    
                    best_match = None
                    best_similarity = 0
                    
                    for result in results:
                        seq = result.get('sequence', {}).get('value', '')
                        if seq:
                            # Check for exact match first
                            if seq == clean_seq:
                                return {
                                    "accession": result.get('primaryAccession', ''),
                                    "protein_name": result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown'),
                                    "organism": result.get('organism', {}).get('scientificName', 'Unknown'),
                                    "similarity": 1.0,
                                    "length": result.get('sequence', {}).get('length', 0)
                                }
                            
                            # Check high similarity
                            if len(seq) == len(clean_seq):
                                matches = sum(1 for i in range(len(seq)) if seq[i] == clean_seq[i])
                                similarity = matches / len(seq)
                                
                                if similarity > best_similarity and similarity > 0.9:
                                    best_similarity = similarity
                                    best_match = result
                    
                    if best_match:
                        return {
                            "accession": best_match.get('primaryAccession', ''),
                            "protein_name": best_match.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown'),
                            "organism": best_match.get('organism', {}).get('scientificName', 'Unknown'),
                            "similarity": best_similarity,
                            "length": best_match.get('sequence', {}).get('length', 0)
                        }
            except Exception as e:
                pass  # Continue to strategy 5
        
        # Strategy 5: Check if this might be a kinase domain from a well-known protein
        if is_likely_kinase:
            # Check against common kinase proteins that might have domain matches
            known_kinases = ['P00533']  # EGFR
            
            for accession in known_kinases:
                try:
                    domain_params = {
                        'query': f'accession:{accession}',
                        'format': 'json',
                        'size': 1,
                        'fields': 'accession,protein_name,organism_name,length,sequence'
                    }
                    
                    response = requests.get(url, params=domain_params, timeout=15)
                    if response.status_code == 200:
                        data = response.json()
                        results = data.get('results', [])
                        
                        if results:
                            result = results[0]
                            full_seq = result.get('sequence', {}).get('value', '')
                            
                            if full_seq:
                                # Check if input sequence is a substring (exact domain match)
                                if clean_seq in full_seq:
                                    start_pos = full_seq.find(clean_seq)
                                    return {
                                        "accession": result.get('primaryAccession', ''),
                                        "protein_name": f"{result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')} (kinase domain)",
                                        "organism": result.get('organism', {}).get('scientificName', 'Unknown'),
                                        "similarity": 1.0,
                                        "length": result.get('sequence', {}).get('length', 0),
                                        "domain_info": f"Kinase domain (residues {start_pos + 1}-{start_pos + len(clean_seq)})",
                                        "domain_start": start_pos + 1
                                    }
                                
                                # Check for high similarity as a potential domain match
                                best_match_score = 0
                                best_match_pos = 0
                                
                                for i in range(len(full_seq) - len(clean_seq) + 1):
                                    subseq = full_seq[i:i+len(clean_seq)]
                                    matches = sum(1 for j in range(len(clean_seq)) if clean_seq[j] == subseq[j])
                                    score = matches / len(clean_seq)
                                    
                                    if score > best_match_score:
                                        best_match_score = score
                                        best_match_pos = i
                                
                                if best_match_score > 0.8:  # High similarity threshold for domain match
                                    return {
                                        "accession": result.get('primaryAccession', ''),
                                        "protein_name": f"{result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Unknown')} (modified kinase domain)",
                                        "organism": result.get('organism', {}).get('scientificName', 'Unknown'),
                                        "similarity": best_match_score,
                                        "length": result.get('sequence', {}).get('length', 0),
                                        "domain_info": f"Similar to kinase domain region (residues {best_match_pos + 1}-{best_match_pos + len(clean_seq)})",
                                        "domain_start": best_match_pos + 1
                                    }
                except Exception as e:
                    continue

        return {"error": "No protein found in UniProt"}

    except Exception as e:
        return {"error": f"Error searching protein: {str(e)}"}


def create_sequence_mapping(input_seq: str, canonical_seq: str, chain_start: int) -> Dict[int, int]:
    """Create mapping between canonical and input sequence positions using sliding window alignment."""
    mapping = {}

    if not canonical_seq:
        # Simple 1:1 mapping if no canonical sequence
        for i in range(len(input_seq)):
            mapping[i + 1] = chain_start + i
        return mapping

    # Find best alignment by sliding input sequence along canonical sequence
    # Try all possible starting positions of input_seq within canonical_seq
    best_start_pos = 0
    best_matches = 0

    # Search range: input sequence could start anywhere in the canonical sequence
    max_search_range = len(canonical_seq) - len(input_seq) + 1

    for start_pos in range(max(0, -len(input_seq) + 1), min(len(canonical_seq), len(canonical_seq) + len(input_seq))):
        matches = 0
        # Count matches when input_seq is positioned at start_pos in canonical_seq
        for i in range(len(input_seq)):
            canonical_idx = start_pos + i
            if 0 <= canonical_idx < len(canonical_seq):
                if canonical_seq[canonical_idx] == input_seq[i]:
                    matches += 1

        if matches > best_matches:
            best_matches = matches
            best_start_pos = start_pos

    # Create mapping based on best alignment
    # Canonical position = best_start_pos + input_position
    for i in range(len(input_seq)):
        canonical_pos = best_start_pos + i + 1  # +1 for 1-based indexing
        if 1 <= canonical_pos <= len(canonical_seq):
            mapping[canonical_pos] = chain_start + i

    return mapping


def map_canonical_to_input_position(canonical_pos: int, mapping: Dict[int, int], chain_start: int) -> Optional[int]:
    """Map canonical position to input sequence position."""
    if canonical_pos in mapping:
        return mapping[canonical_pos]

    # Fallback: assume 1:1 mapping if not found
    return canonical_pos - 1 + chain_start if canonical_pos > 0 else None


def extract_mutation_residues(description: str, canonical_seq: str, position: int) -> Tuple[Optional[str], Optional[str]]:
    """Extract wild-type and mutant amino acids from variant description."""
    # Try different patterns
    patterns = [
        r'([A-Z])\s*->\s*([A-Z])',  # A -> V
        r'([A-Z])(\d+)([A-Z])',      # A50V
        r'p\.([A-Z])[a-z]*(\d+)([A-Z])[a-z]*',  # p.Ala50Val
        r'([A-Z])[a-z]*\s*->\s*([A-Z])[a-z]*'   # Ala -> Val
    ]

    for pattern in patterns:
        match = re.search(pattern, description.upper())
        if match:
            if len(match.groups()) >= 3:  # A50V pattern
                return match.group(1), match.group(3)
            else:  # A -> V pattern
                return match.group(1), match.group(2)

    # Fallback: try to extract from canonical sequence
    if canonical_seq and 1 <= position <= len(canonical_seq):
        wt_residue = canonical_seq[position - 1]
        # Try to find mutant residue in description
        mut_match = re.search(r'([A-Z])', description.replace(wt_residue, '', 1))
        if mut_match:
            return wt_residue, mut_match.group(1)

    return None, None


def categorize_amino_acid_change(wt_residue: str, mut_residue: str) -> str:
    """Categorize the type of amino acid change."""
    if not wt_residue or not mut_residue:
        return "Unknown"
    
    # Amino acid properties
    hydrophobic = set('AILMFPWV')
    polar = set('NQSTY')
    positive = set('KRH')
    negative = set('DE')
    aromatic = set('FWY')
    small = set('AGS')
    
    wt = wt_residue.upper()
    mut = mut_residue.upper()
    
    if wt == mut:
        return "Silent"
    
    # Property-based categorization
    if wt in positive and mut in negative:
        return "Charge reversal"
    elif wt in negative and mut in positive:
        return "Charge reversal"
    elif wt in hydrophobic and mut in polar:
        return "Hydrophobic to polar"
    elif wt in polar and mut in hydrophobic:
        return "Polar to hydrophobic"
    elif wt in aromatic and mut not in aromatic:
        return "Loss of aromaticity"
    elif wt not in aromatic and mut in aromatic:
        return "Gain of aromaticity"
    elif wt in small and mut not in small:
        return "Size increase"
    elif wt not in small and mut in small:
        return "Size decrease"
    else:
        return "Conservative"


def predict_functional_severity(wt_residue: str, mut_residue: str, change_type: str) -> str:
    """Predict the functional severity of a mutation."""
    if change_type == "Silent":
        return "None"
    elif change_type in ["Charge reversal"]:
        return "High"
    elif change_type in ["Hydrophobic to polar", "Polar to hydrophobic", "Loss of aromaticity"]:
        return "Moderate"
    elif change_type in ["Size increase", "Size decrease"]:
        return "Low"
    else:
        return "Low"


def predict_structural_impact(change_type: str, severity: str) -> str:
    """Predict structural impact based on change type and severity."""
    if change_type == "Silent":
        return "None"
    elif severity == "High":
        return "Likely disrupts local structure"
    elif severity == "Moderate":
        return "May alter local conformation"
    else:
        return "Minimal structural change"


def query_mutations_for_protein_detailed(protein_info: Dict, input_sequence: str, chain_start: int = 1, max_mutations: int = 20) -> List[Dict]:
    """
    Query known mutations for a protein from public databases with detailed information.
    Returns a list of mutation dictionaries with comprehensive details.
    """
    mutations = []

    if "error" in protein_info:
        return []

    try:
        accession = protein_info.get("accession", "")
        if not accession:
            return []

        # Clean input sequence for position mapping
        clean_input_seq = re.sub(r'[^A-Z]', '', input_sequence.upper())
        
        # Check if this is a domain sequence with known start position
        domain_start = protein_info.get("domain_start", None)

        # Query UniProt for variants of this protein
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f'accession:{accession}',
            'format': 'json',
            'size': 1,
            'fields': 'accession,sequence,ft_variant'
        }

        response = requests.get(url, params=params, timeout=15)
        if response.status_code == 200:
            data = response.json()
            results = data.get('results', [])

            if results:
                result = results[0]
                canonical_seq = result.get('sequence', {}).get('value', '')
                features = result.get('features', [])

                # Handle domain sequences vs full sequences differently
                if domain_start:
                    # For domain sequences, map positions directly to the domain region
                    domain_end = domain_start + len(clean_input_seq) - 1
                    
                    for feature in features:
                        if feature.get('type') == 'Natural variant':
                            # Extract mutation details
                            location = feature.get('location', {})
                            canonical_pos = location.get('start', {}).get('value')
                            description = feature.get('description', '')
                            evidence = feature.get('evidences', [])

                            if canonical_pos and description:
                                # Check if mutation is within the domain region
                                if domain_start <= canonical_pos <= domain_end:
                                    # Calculate position within the input domain sequence
                                    input_pos = canonical_pos - domain_start + 1

                                    # Extract wild-type and mutant residues
                                    wt_residue, mut_residue = extract_mutation_residues(description, canonical_seq, canonical_pos)

                                    if wt_residue and mut_residue:
                                        # For domain sequences, we expect the canonical sequence to match
                                        # at the domain position, not the input position
                                        if 1 <= canonical_pos <= len(canonical_seq):
                                            canonical_residue = canonical_seq[canonical_pos - 1]
                                            if canonical_residue == wt_residue:
                                                # Create mutation in input sequence numbering
                                                input_mutation_pos = chain_start + input_pos - 1

                                                # Categorize the change
                                                change_type = categorize_amino_acid_change(wt_residue, mut_residue)
                                                severity = predict_functional_severity(wt_residue, mut_residue, change_type)
                                                structural_impact = predict_structural_impact(change_type, severity)

                                                mutation_data = {
                                                    'mutation_code': f"{wt_residue}{input_mutation_pos}{mut_residue}",
                                                    'canonical_code': f"{wt_residue}{canonical_pos}{mut_residue}",
                                                    'input_position': input_mutation_pos,
                                                    'canonical_position': canonical_pos,
                                                    'wt_residue': wt_residue,
                                                    'mut_residue': mut_residue,
                                                    'change_type': change_type,
                                                    'severity': severity,
                                                    'structural_impact': structural_impact,
                                                    'evidence_count': len(evidence),
                                                    'description': description,
                                                    'uniprot_accession': accession
                                                }

                                                mutations.append(mutation_data)

                else:
                    # Original logic for full sequences
                    # Create sequence alignment mapping
                    seq_mapping = create_sequence_mapping(clean_input_seq, canonical_seq, chain_start)

                    for feature in features:
                        if feature.get('type') == 'Natural variant':
                            # Extract mutation details
                            location = feature.get('location', {})
                            canonical_pos = location.get('start', {}).get('value')
                            description = feature.get('description', '')
                            evidence = feature.get('evidences', [])

                            if canonical_pos and description:
                                # Map canonical position to input sequence position
                                input_pos = map_canonical_to_input_position(canonical_pos, seq_mapping, chain_start)

                                if input_pos and 1 <= input_pos <= len(clean_input_seq):
                                    # Extract wild-type and mutant residues
                                    wt_residue, mut_residue = extract_mutation_residues(description, canonical_seq, canonical_pos)

                                    if wt_residue and mut_residue:
                                        # Verify wild-type residue matches input sequence
                                        if clean_input_seq[input_pos - 1] == wt_residue:
                                            # Categorize the change
                                            change_type = categorize_amino_acid_change(wt_residue, mut_residue)
                                            severity = predict_functional_severity(wt_residue, mut_residue, change_type)
                                            structural_impact = predict_structural_impact(change_type, severity)

                                            mutation_data = {
                                                'mutation_code': f"{wt_residue}{input_pos}{mut_residue}",
                                                'canonical_code': f"{wt_residue}{canonical_pos}{mut_residue}",
                                                'input_position': input_pos,
                                                'canonical_position': canonical_pos,
                                                'wt_residue': wt_residue,
                                                'mut_residue': mut_residue,
                                                'change_type': change_type,
                                                'severity': severity,
                                                'structural_impact': structural_impact,
                                                'evidence_count': len(evidence),
                                                'description': description,
                                                'uniprot_accession': accession
                                            }

                                            mutations.append(mutation_data)

    except Exception as e:
        st.error(f"Error querying mutations: {e}")

    return mutations


def format_mutations_for_input(mutations: List[Dict]) -> str:
    """Format mutations as a simple comma-separated list for user input."""
    if not mutations:
        return ""

    mutation_codes = [m['mutation_code'] for m in mutations if 'mutation_code' in m]
    return ", ".join(mutation_codes[:10])  # Limit to first 10 mutations


def format_mutations_table_improved(mutations: List[Dict]) -> str:
    """Format mutations as a readable text table without redundant position columns."""
    if not mutations:
        return "No mutations found."

    # Calculate maximum widths for each column
    mutation_width = max(len("Mutation"), max(len(m['mutation_code']) for m in mutations)) + 1
    change_width = max(len("Change Type"), max(len(m.get('change_type', '')) for m in mutations)) + 1
    severity_width = max(len("Severity"), max(len(m.get('severity', '')) for m in mutations)) + 1
    impact_width = max(len("Structural Impact"), max(len(m.get('structural_impact', '')) for m in mutations)) + 1
    evidence_width = max(len("Evidence"), max(len(str(m.get('evidence_count', 0))) for m in mutations)) + 1

    # Create header
    header = (f"{'Mutation':<{mutation_width}} "
              f"{'Change Type':<{change_width}} "
              f"{'Severity':<{severity_width}} "
              f"{'Structural Impact':<{impact_width}} "
              f"{'Evidence':<{evidence_width}}")
    
    separator = "-" * len(header)
    
    # Create rows
    rows = []
    for m in mutations:
        row = (f"{m['mutation_code']:<{mutation_width}} "
               f"{m.get('change_type', ''):<{change_width}} "
               f"{m.get('severity', ''):<{severity_width}} "
               f"{m.get('structural_impact', ''):<{impact_width}} "
               f"{m.get('evidence_count', 0):<{evidence_width}}")
        rows.append(row)
    
    return f"{header}\n{separator}\n" + "\n".join(rows)


def display_mutation_discovery_section(wt_protein_sequence: str = None, chain_start: int = 1, max_mutations: int = 50) -> Optional[str]:
    """
    Display mutation discovery section and return formatted mutation string.

    Args:
        wt_protein_sequence: Wild-type protein sequence
        chain_start: Starting residue number for the chain
        max_mutations: Maximum number of mutations to discover (default: 50)

    Returns:
        Formatted mutation string for input or None if not discovered
    """
    
    # Get latest protein sequence
    latest_protein_seq = None
    if 'wt_protein_input' in st.session_state and st.session_state['wt_protein_input']:
        latest_protein_seq = st.session_state['wt_protein_input']
    elif wt_protein_sequence:
        latest_protein_seq = wt_protein_sequence
    else:
        latest_protein_seq = ''

    if not latest_protein_seq or not latest_protein_seq.strip():
        st.info("Please set the input mode to 'Mutation Mode' and enter a wild-type protein sequence above to discover mutations.")
        return None

    # Validate sequence
    is_valid, error_msg, chains_dict, clean_seq = validate_protein_sequence(latest_protein_seq)

    if not is_valid:
        st.error(f"Invalid protein sequence: {error_msg}")
        return None

    # Use first chain for mutation discovery
    if chains_dict and len(chains_dict) > 0:
        first_chain_id = sorted(chains_dict.keys())[0]
        first_chain_seq = chains_dict[first_chain_id]
    else:
        first_chain_seq = clean_seq

    # Description text
    st.markdown("Automatically retrieve known variants for the current wild-type protein sequence.")

    # Use unlimited mutations - find as many as possible
    max_mutations_input = 999999  # Effectively unlimited

    # Chain selector - always show dropdown for consistency
    if chains_dict and len(chains_dict) > 0:
        available_chains = sorted(chains_dict.keys())
        selected_chain_for_discovery = st.selectbox(
            "Chain to discover mutations for",
            options=available_chains,
            index=0,
            help="Select which chain to find mutations for",
            key="mutation_chain_selector"
        )
        # Update first_chain_seq based on selection
        first_chain_seq = chains_dict[selected_chain_for_discovery]
        first_chain_id = selected_chain_for_discovery
    else:
        # Single chain - still show dropdown with just one option for consistency
        selected_chain_for_discovery = st.selectbox(
            "Chain to discover mutations for",
            options=['A'],
            index=0,
            help="Select which chain to find mutations for",
            key="mutation_chain_selector"
        )
        first_chain_id = 'A'

    # Centered button
    col_left, col_btn, col_right = st.columns([1, 1, 1])

    with col_btn:
        discover_button = st.button(
            "Discover Mutations",
            key="discover_mutations_btn",
            icon=":material/search:",
            type="secondary",
            use_container_width=True
        )

    # Initialize session state for discovered mutations
    if 'discovered_mutations' not in st.session_state:
        st.session_state.discovered_mutations = None
    if 'discovery_sequence_hash' not in st.session_state:
        st.session_state.discovery_sequence_hash = None

    # Create a hash of the current sequence to detect changes
    import hashlib
    current_seq_hash = hashlib.md5(first_chain_seq.encode()).hexdigest()

    # Clear stored mutations if sequence changed
    if st.session_state.discovery_sequence_hash != current_seq_hash:
        st.session_state.discovered_mutations = None
        st.session_state.discovery_sequence_hash = None

    if discover_button:
        detailed_mutations: List[Dict] = []
        try:
            with st.spinner("Searching protein databases and discovering mutations..."):
                detailed_mutations = discover_mutations_for_sequence(
                    first_chain_seq,
                    max_mutations=max_mutations_input,
                    chain_start=chain_start
                )

            if not detailed_mutations:
                st.warning("No mutations found in public databases for this protein sequence.")
                st.session_state.discovered_mutations = None
                return None

            st.success(f"Found {len(detailed_mutations)} mutations from database search")
            # Store in session state
            st.session_state.discovered_mutations = detailed_mutations
            st.session_state.discovery_sequence_hash = current_seq_hash

        except Exception as e:
            st.info(f"Database search failed, trying standard search: {e}")

            with st.spinner("Searching UniProt for protein information..."):
                protein_info = get_protein_info(first_chain_seq)

            if "error" in protein_info:
                st.error(f"Could not find protein in UniProt: {protein_info['error']}")
                st.session_state.discovered_mutations = None
                return None

            similarity = protein_info.get('similarity', 1.0)
            if similarity >= 1.0:
                st.success(f"Found protein: {protein_info.get('protein_name', 'Unknown')} ({protein_info.get('accession', 'Unknown')})")
            else:
                st.success(f"Found protein: {protein_info.get('protein_name', 'Unknown')} ({protein_info.get('accession', 'Unknown')}) - {similarity:.1%} similarity")

            with st.spinner("Discovering known mutations..."):
                detailed_mutations = query_mutations_for_protein_detailed(
                    protein_info, first_chain_seq, chain_start, max_mutations_input
                )

            if not detailed_mutations:
                st.warning("No mutations found in UniProt for this protein sequence.")
                st.session_state.discovered_mutations = None
                return None

            st.success(f"Found {len(detailed_mutations)} mutations")
            # Store in session state
            st.session_state.discovered_mutations = detailed_mutations
            st.session_state.discovery_sequence_hash = current_seq_hash

    # Display mutations if they exist in session state (persists across reruns)
    if st.session_state.discovered_mutations:
        detailed_mutations = st.session_state.discovered_mutations

        # Filter to only standard substitution mutations and sort by position
        def is_substitution_mutation(mutation_code):
            """Check if mutation is a standard single amino acid substitution (e.g., R127Q)"""
            # Match pattern: Single letter + numbers + Single letter
            import re
            return bool(re.match(r'^[A-Z]\d+[A-Z]$', mutation_code))

        def get_mutation_position(mutation_code):
            """Extract numeric position from mutation code (e.g., R127Q -> 127)"""
            import re
            match = re.search(r'\d+', mutation_code)
            return int(match.group()) if match else 0

        # Filter and sort mutations
        filtered_mutations = [m for m in detailed_mutations if is_substitution_mutation(m['mutation_code'])]
        filtered_mutations.sort(key=lambda x: get_mutation_position(x['mutation_code']))

        if not filtered_mutations:
            st.warning("No standard substitution mutations found")
            return None


        # Initialize session state for filters
        if 'mutation_offset' not in st.session_state:
            st.session_state.mutation_offset = 0
        if 'mutation_range_start' not in st.session_state:
            st.session_state.mutation_range_start = None
        if 'mutation_range_end' not in st.session_state:
            st.session_state.mutation_range_end = None

        # Calculate suggested offset
        sample_mutation = filtered_mutations[0] if filtered_mutations else None
        if sample_mutation:
            canonical_start = sample_mutation.get('canonical_position', 0)
            user_pos = sample_mutation.get('position_in_user_sequence', 1)
            suggested_offset = canonical_start - user_pos if canonical_start and user_pos else 0
        else:
            suggested_offset = 0

        # Apply filters and offset
        display_mutations = []
        for m in filtered_mutations:
            # Extract position from mutation_code
            import re
            match = re.search(r'\d+', m['mutation_code'])
            if match:
                original_pos = int(match.group())
                # Apply offset
                adjusted_pos = original_pos + st.session_state.mutation_offset

                # Apply range filter
                if st.session_state.mutation_range_start and adjusted_pos < st.session_state.mutation_range_start:
                    continue
                if st.session_state.mutation_range_end and adjusted_pos > st.session_state.mutation_range_end:
                    continue

                # Create adjusted mutation with new numbering
                adjusted_mutation = m.copy()
                wt = m['mutation_code'][0]
                mut = m['mutation_code'][-1]
                adjusted_mutation['mutation_code'] = f"{wt}{adjusted_pos}{mut}"
                display_mutations.append(adjusted_mutation)

        # Show info about current settings
        if st.session_state.mutation_offset != 0 or st.session_state.mutation_range_start or st.session_state.mutation_range_end:
            filter_info = []
            if st.session_state.mutation_offset != 0:
                filter_info.append(f"Offset: +{st.session_state.mutation_offset}")
            if st.session_state.mutation_range_start:
                filter_info.append(f"From position {st.session_state.mutation_range_start}")
            if st.session_state.mutation_range_end:
                filter_info.append(f"To position {st.session_state.mutation_range_end}")
            st.caption(f"Active filters: {', '.join(filter_info)} | Showing {len(display_mutations)}/{len(filtered_mutations)} mutations")

        filtered_mutations = display_mutations if display_mutations else filtered_mutations

        # Filter controls section - moved above the table
        st.markdown(":material/tune: **Adjust Numbering & Filters**")

        # First row: inputs
        filter_col1, filter_col2, filter_col3, btn_col1, btn_col2 = st.columns([1, 1, 1, 0.6, 0.6])

        with filter_col1:
            offset_input = st.number_input(
                "Position Offset",
                value=suggested_offset,
                help=f"Adjust residue numbering. Suggested: {suggested_offset} (based on UniProt alignment)",
                key="mutation_offset_input"
            )

        with filter_col2:
            range_start = st.number_input(
                "Start Position",
                value=None,
                min_value=1,
                help="Only show mutations starting from this position (optional)",
                key="mutation_range_start_input"
            )

        with filter_col3:
            range_end = st.number_input(
                "End Position",
                value=None,
                min_value=1,
                help="Only show mutations up to this position (optional)",
                key="mutation_range_end_input"
            )

        # Buttons inline with End Position
        with btn_col1:
            if st.button(":material/check: Apply", key="apply_mutation_filters", type="tertiary", use_container_width=True):
                st.session_state.mutation_offset = offset_input
                st.session_state.mutation_range_start = range_start
                st.session_state.mutation_range_end = range_end
                st.rerun()
        with btn_col2:
            if st.button(":material/refresh: Reset", key="reset_mutation_filters", type="tertiary", use_container_width=True):
                st.session_state.mutation_offset = suggested_offset
                st.session_state.mutation_range_start = None
                st.session_state.mutation_range_end = None
                st.rerun()

        # Display mutation analysis table (full width)
        st.markdown("##### :material/analytics: Mutation Analysis Table")

        # Create DataFrame with prediction scores
        mutation_df = pd.DataFrame([
            {
                "Mutation (Your Seq)": m['mutation_code'],
                "Canonical Position": m.get('canonical_code', m['mutation_code']),
                "SIFT Score": m.get('sift_score', '') if m.get('sift_score') is not None else '-',
                "SIFT Prediction": m.get('sift_prediction', '-'),
                "PolyPhen Score": m.get('polyphen_score', '') if m.get('polyphen_score') is not None else '-',
                "PolyPhen Prediction": m.get('polyphen_prediction', '-'),
                "Clinical Significance": m.get('clinical_significance', '-'),
                "References": m.get('xrefs', [])
            }
            for m in filtered_mutations
        ])

        # Format references (renamed from literature)
        def format_references(xrefs):
            if not xrefs or len(xrefs) == 0:
                return '-'
            # Extract PMIDs or other references
            refs = []
            for xref in xrefs[:3]:  # Show max 3 references
                if isinstance(xref, dict):
                    db = xref.get('database', '')
                    ref_id = xref.get('id', '')
                    if 'pubmed' in db.lower() or 'pmid' in db.lower():
                        refs.append(f"PMID:{ref_id}")
                    else:
                        refs.append(f"{db}:{ref_id}")
                elif isinstance(xref, str):
                    refs.append(xref)
            return ', '.join(refs) if refs else '-'

        mutation_df['References'] = mutation_df['References'].apply(format_references)

        st.dataframe(
            mutation_df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Mutation (Your Seq)": st.column_config.TextColumn("Mutation (Your Seq)", width="small", help="Mutation using your sequence numbering (adjusted for chain start position)"),
                "Canonical Position": st.column_config.TextColumn("Canonical", width="small", help="Mutation in canonical UniProt numbering"),
                "SIFT Score": st.column_config.TextColumn("SIFT Score", width="small"),
                "SIFT Prediction": st.column_config.TextColumn("SIFT", width="small"),
                "PolyPhen Score": st.column_config.TextColumn("PolyPhen Score", width="small"),
                "PolyPhen Prediction": st.column_config.TextColumn("PolyPhen", width="small"),
                "Clinical Significance": st.column_config.TextColumn("Clinical Significance", width="medium"),
                "References": st.column_config.TextColumn("References", width="medium")
            }
        )

        # Format mutations for input (display below table)
        st.markdown("##### :material/content_copy: Copy-Pastable Mutation List")
        formatted_mutations = format_mutations_for_input(filtered_mutations)
        st.code(formatted_mutations, language=None)

        # Column explanations
        st.markdown("##### :material/info: Column Descriptions")
        st.markdown("""
        - **Mutation (Your Seq)**: Mutation using your sequence numbering (with applied offset). **Use this for mutations input.**
        - **Canonical**: Mutation in standard UniProt numbering - for reference lookup
        - **SIFT**: Predicts tolerance - deleterious if score ≤0.05 (may be empty if not available in database)
        - **PolyPhen**: Predicts damage - probably/possibly damaging or benign (may be empty if not available)
        - **Clinical Significance**: Known clinical impact from databases
        - **References**: References (PMID or database IDs) supporting this variant

        **Numbering Adjustment:** Use the "Adjust Numbering & Filters" controls above to set a position offset if your sequence
        numbering differs from the canonical UniProt positions. The suggested offset aligns your sequence with the full protein.
        Apply filters to show only mutations in a specific residue range.
        """)

        return formatted_mutations

    return None
