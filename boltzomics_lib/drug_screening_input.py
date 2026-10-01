from typing import Any, Dict, List, Optional, Tuple
import streamlit as st
import re
import math
from functools import lru_cache
import pandas as pd
import requests
import time

from mutation_discovery import discover_mutations_for_sequence

try:
    from alphamissense_lookup import (
        annotate_mutations_with_alphamissense,
        get_default_alphamissense_lookup,
    )
except ImportError:
    annotate_mutations_with_alphamissense = None
    get_default_alphamissense_lookup = None

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

def _deduplicate_name(name: str, used: set, prefix: str) -> str:
    base = re.sub(r"\s+", " ", str(name or "").strip()) or prefix
    candidate = base
    suffix = 2
    while candidate.casefold() in used:
        candidate = f"{base}_{suffix}"
        suffix += 1
    used.add(candidate.casefold())
    return candidate


def _normalize_sequence_fragment(fragment: str) -> Tuple[str, List[str]]:
    corrections: List[str] = []
    original = fragment
    fragment = fragment.replace("\ufeff", "")
    if re.search(r"\d", fragment):
        fragment = re.sub(r"\d", "", fragment)
        corrections.append("removed pasted residue numbers")
    compact = re.sub(r"\s+", "", fragment)
    if compact != fragment or "\n" in original or "\r" in original:
        corrections.append("removed sequence whitespace/line breaks")
    if compact.endswith("*") and "*" not in compact[:-1]:
        compact = compact[:-1]
        corrections.append("removed terminal FASTA stop marker (*)")
    return compact.upper(), corrections


def parse_fasta_sequences(fasta_text: str, return_corrections: bool = False):
    """
    Parse FASTA format text and return list of (name, sequence) tuples.
    
    Args:
        fasta_text: FASTA format text with sequences
        
    Returns:
        List of tuples containing (sequence_name, sequence)
    """
    sequences: List[Tuple[str, str]] = []
    corrections: List[str] = []
    used_names: set = set()
    current_name = None
    current_fragments: List[str] = []
    
    lines = str(fasta_text or "").replace("\r\n", "\n").replace("\r", "\n").split("\n")

    def flush_sequence() -> None:
        nonlocal current_name, current_fragments
        if not current_fragments:
            return
        sequence, applied = _normalize_sequence_fragment("\n".join(current_fragments))
        if sequence:
            name = _deduplicate_name(current_name, used_names, f"Protein_{len(sequences) + 1}")
            if current_name is None:
                corrections.append("added a missing FASTA identifier")
            elif name != str(current_name).strip():
                corrections.append("made duplicate/blank FASTA identifiers unique")
            sequences.append((name, sequence))
            corrections.extend(applied)
        current_fragments = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(";"):
            corrections.append("ignored FASTA comment lines")
            continue
            
        if line.startswith('>'):
            flush_sequence()
            current_name = line[1:].strip() or None
        else:
            current_fragments.append(line)
    
    # Add the last sequence
    flush_sequence()
    unique_corrections = list(dict.fromkeys(corrections))
    return (sequences, unique_corrections) if return_corrections else sequences

def parse_smiles_list(smiles_text: str, return_corrections: bool = False):
    """
    Parse SMILES text in FASTA format and return list of (name, smiles) tuples.
    
    Args:
        smiles_text: FASTA format text with SMILES strings
        
    Returns:
        List of tuples containing (drug_name, smiles)
    """
    drugs: List[Tuple[str, str]] = []
    corrections: List[str] = []
    used_names: set = set()
    current_name = None
    current_fragments: List[str] = []

    lines = str(smiles_text or "").replace("\r\n", "\n").replace("\r", "\n").split("\n")
    has_headers = any(line.strip().startswith(">") for line in lines)

    def clean_smiles(value: str) -> str:
        cleaned = value.strip().strip("`\"'")
        compact = re.sub(r"\s+", "", cleaned)
        if compact != value.strip():
            corrections.append("removed SMILES whitespace or wrapping quotes")
        return compact

    if not has_headers:
        for raw_line in lines:
            line = raw_line.strip()
            if not line:
                continue
            name = None
            value = line
            for delimiter in ("\t", ","):
                if delimiter in line:
                    possible_name, possible_value = line.split(delimiter, 1)
                    if possible_name.strip() and possible_value.strip():
                        name, value = possible_name.strip(), possible_value.strip()
                        corrections.append("converted a name/SMILES row to FASTA-style input")
                        break
            name = _deduplicate_name(name, used_names, f"Drug_{len(drugs) + 1}")
            if name.startswith("Drug_"):
                corrections.append("added missing drug identifiers")
            drugs.append((name, clean_smiles(value)))
        unique_corrections = list(dict.fromkeys(corrections))
        return (drugs, unique_corrections) if return_corrections else drugs

    def flush_drug() -> None:
        nonlocal current_name, current_fragments
        if current_fragments:
            value = clean_smiles("".join(current_fragments))
            if value:
                name = _deduplicate_name(current_name, used_names, f"Drug_{len(drugs) + 1}")
                if current_name is None:
                    corrections.append("added a missing drug identifier")
                drugs.append((name, value))
        current_fragments = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            flush_drug()
            current_name = line[1:].strip() or None
        else:
            current_fragments.append(line)
    
    # Add the last drug
    flush_drug()
    unique_corrections = list(dict.fromkeys(corrections))
    return (drugs, unique_corrections) if return_corrections else drugs

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
    raw_input = str(protein_seq or "")
    if any(line.strip().startswith(">") for line in raw_input.splitlines()):
        parsed = parse_fasta_sequences(raw_input)
        if len(parsed) != 1:
            return False, "Expected exactly one FASTA record in this field.", {}, raw_input
        protein_seq = parsed[0][1]
    else:
        protein_seq, _ = _normalize_sequence_fragment(raw_input)
    if not protein_seq.strip():
        return False, "Protein sequence is empty after normalization.", {}, protein_seq
    
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

def parse_pocket_residue_positions(chain: str, positions: str, chain_length: Optional[int] = None) -> List[List[Any]]:
    """Expand friendly position lists and ranges, using input-sequence numbering."""
    residues = []
    positions = re.sub(r"\s*[-–—]\s*", "-", positions.strip())
    for part in re.split(r"[\s,]+", positions):
        if not part:
            continue
        match = re.fullmatch(r"([1-9]\d*)(?:-([1-9]\d*))?", part)
        if not match:
            raise ValueError(f"Invalid residue position '{part}'. Use numbers or ranges, for example 123, 125, 130-135.")
        start = int(match[1])
        end = int(match[2] or start)
        if end < start:
            raise ValueError(f"Residue range '{part}' must go from a smaller to a larger position.")
        if chain_length is not None and end > chain_length:
            raise ValueError(f"Chain {chain} has {chain_length} residues; position {end} is outside the input sequence.")
        if end - start > 10000:
            raise ValueError("Select a pocket range of at most 10001 residues.")
        residues.extend(range(start, end + 1))
    return [[chain, residue] for residue in dict.fromkeys(residues)]


def normalize_binding_pocket_constraints(constraints: Any) -> List[Dict[str, Any]]:
    """Accept legacy single pockets and multiple independent Boltz2 pockets."""
    if constraints is None:
        return []
    if isinstance(constraints, dict):
        constraints = [constraints] if constraints.get("contacts") else []
    if not isinstance(constraints, list):
        raise ValueError("Binding pocket constraints must be an object or a list of objects.")
    normalized = []
    for index, pocket in enumerate(constraints, 1):
        if not isinstance(pocket, dict):
            raise ValueError(f"Pocket constraint {index} must be an object.")
        binder = str(pocket.get("binder", "X")).strip()
        if not binder:
            raise ValueError(f"Pocket constraint {index} needs a binder chain.")
        contacts = pocket.get("contacts")
        if not isinstance(contacts, list) or not contacts:
            raise ValueError(f"Pocket constraint {index} needs at least one pocket residue or ligand atom.")
        # Contact and pocket endpoints use the same Boltz token syntax.
        endpoints = normalize_distance_constraints([
            {"token1": contact, "token2": contact,
             "max_distance": pocket.get("max_distance", 6.0),
             "force": pocket.get("force", False)}
            for contact in contacts
        ])
        normalized.append({
            **pocket,
            "binder": binder,
            "contacts": [list(contact) for contact in dict.fromkeys(
                tuple(endpoint["token1"]) for endpoint in endpoints
            )],
            "max_distance": endpoints[0]["max_distance"],
            "force": endpoints[0]["force"],
        })
    seen_contacts = {}
    for index, pocket in enumerate(normalized, 1):
        for chain, target in pocket["contacts"]:
            key = (pocket["binder"], chain, target)
            if key in seen_contacts and seen_contacts[key] != pocket["max_distance"]:
                raise ValueError(f"Pocket constraint {index}: {chain},{target} already has a different distance limit for binder {pocket['binder']}. Use one limit for this contact.")
            seen_contacts[key] = pocket["max_distance"]
    return normalized


def binding_pocket_yaml_entries(constraints: Any) -> List[Dict[str, Any]]:
    """Keep editor metadata out of the Boltz YAML schema."""
    return [{"pocket": {key: pocket[key] for key in
                         ("binder", "contacts", "max_distance", "force")}}
            for pocket in normalize_binding_pocket_constraints(constraints)]


def display_applied_pocket_constraints(constraints: Any) -> None:
    pockets = normalize_binding_pocket_constraints(constraints)
    if pockets:
        st.markdown("**Applied ligand → pocket constraints**")
        for index, pocket in enumerate(pockets, 1):
            contacts = ", ".join(f"{chain},{residue}" for chain, residue in pocket["contacts"])
            enforcement = "potential enabled" if pocket["force"] else "model guidance"
            st.write(
                f"{index}. Binder {pocket['binder']} → each of [{contacts}] "
                f"≤ {pocket['max_distance']:g} Å · {enforcement}"
            )


@lru_cache(maxsize=256)
def boltz_ligand_atoms(smiles: str, affinity: bool = True) -> Tuple[Tuple[str, str], ...]:
    """Use Boltz 2.2.1's SMILES preparation and canonical atom naming rules."""
    from rdkit import Chem
    if affinity:
        from chembl_structure_pipeline.exclude_flag import exclude_flag
        from chembl_structure_pipeline.standardizer import standardize_mol
        from rdkit.Chem.MolStandardize import rdMolStandardize
        molecule = Chem.MolFromSmiles(smiles, sanitize=False)
        if molecule is None or exclude_flag(molecule, includeRDKitSanitization=False):
            raise ValueError("This ligand cannot be prepared for Boltz affinity prediction.")
        molecule = rdMolStandardize.LargestFragmentChooser().choose(molecule)
        smiles = Chem.MolToSmiles(standardize_mol(molecule))
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("The ligand SMILES could not be parsed.")
    molecule = Chem.AddHs(molecule)
    ranks = Chem.CanonicalRankAtoms(molecule)
    return tuple((atom.GetSymbol().upper() + str(rank + 1), atom.GetSymbol())
                 for atom, rank in zip(molecule.GetAtoms(), ranks) if atom.GetAtomicNum() != 1)


def selected_residue_rows(contacts: List[List[Any]], chains: Dict[str, str], limit: float) -> List[Dict[str, Any]]:
    """Build a deduplicated confirmation table from input-sequence positions."""
    rows = []
    seen = set()
    for chain, position in contacts:
        try:
            position = int(position)
        except (TypeError, ValueError):
            continue
        if (chain, position) in seen:
            continue
        seen.add((chain, position))
        sequence = chains.get(chain, "")
        residue = sequence[position - 1] if 1 <= position <= len(sequence) else "Unavailable"
        rows.append({"Chain": chain, "Position (1-based)": position,
                     "Residue": residue, "Ligand distance limit (Å)": limit})
    return rows


def validate_boltz_constraint_targets(constraints: List[Dict[str, Any]], sequences: List[Dict[str, Any]]) -> None:
    """Check generated constraints against the chains actually written to the YAML."""
    chains = {}
    for entry in sequences:
        kind, data = next(iter(entry.items()))
        ids = data["id"] if isinstance(data["id"], list) else [data["id"]]
        for chain in ids:
            if chain in chains:
                raise ValueError(f"Duplicate chain ID {chain} in the prediction input.")
            chains[chain] = (kind, data)

    def validate_endpoint(token):
        chain, target = token
        if chain not in chains:
            raise ValueError(f"Constraint references chain {chain}, which is missing from this prediction.")
        kind, data = chains[chain]
        if kind == "ligand":
            if not isinstance(target, str):
                raise ValueError(f"Ligand chain {chain} needs an atom name, not a residue number.")
            if data.get("smiles"):
                names = {name for name, _ in boltz_ligand_atoms(data["smiles"], affinity=chain == "X")}
                if target not in names:
                    raise ValueError(f"Ligand atom {chain},{target} does not exist in this molecule. Available atoms: {', '.join(sorted(names))}.")
        elif not isinstance(target, int) or not 1 <= target <= len(data["sequence"]):
            raise ValueError(f"Chain {chain} needs a residue position between 1 and {len(data['sequence'])}.")

    for constraint in constraints:
        if "pocket" in constraint:
            pocket = constraint["pocket"]
            if pocket["binder"] not in chains:
                raise ValueError(f"Pocket binder {pocket['binder']} is missing from this prediction. A ligand pocket requires a ligand input.")
            for contact in pocket["contacts"]:
                if contact[0] == pocket["binder"]:
                    raise ValueError("Pocket contacts must belong to a different chain from the binder.")
                validate_endpoint(contact)
        if "contact" in constraint:
            contact = constraint["contact"]
            if contact["token1"] == contact["token2"]:
                raise ValueError("A pair constraint needs two different endpoints.")
            validate_endpoint(contact["token1"])
            validate_endpoint(contact["token2"])


def _remove_editor_row(rows_key: str, row_id: int) -> None:
    # Callbacks run before rendering, so surviving widgets keep their values.
    st.session_state[rows_key] = [row for row in st.session_state[rows_key] if row["id"] != row_id]


def _add_pocket_chain(groups_key: str) -> None:
    next_id_key = f"{groups_key}_next_id"
    next_id = st.session_state[next_id_key]
    st.session_state[groups_key].append({"id": next_id, "chain": "A", "positions": ""})
    st.session_state[next_id_key] = next_id + 1


def display_constraint_examples() -> None:
    st.markdown("Choose what you want to bring closer: **the ligand to a pocket**, or **two specific endpoints**.")
    with st.expander("Examples and how to use", expanded=False):
        st.markdown("""
**Place the ligand near a pocket:** Under Ligand → binding pocket, select protein chain A,
enter residue positions `123, 125, 130-135`, and set the limit to **6 Å**.
This guides the ligand within 6 Å of each selected residue. Click **Apply pocket constraints**.

**Bring two protein residues closer:** Under Pair distance, choose Protein residue for both endpoints.
Select chain A, position 123 and chain B, position 45. Set **6 Å** and click **Apply pair constraints**.

**Target a particular ligand atom:** Choose Protein residue A,123 and Ligand atom X,O1.
The limit applies between that residue and that named ligand atom.
Replace O1 with the exact atom name from your ligand's Boltz structure.

**Use different limits:** Add one pocket constraint for A,123 at **6 Å** and another for A,125 at **8 Å**.
You can also add pair constraints. All applied constraints are combined in the prediction;
multiple pockets are simultaneous requirements, not alternative binding sites.

**Give each pair its own limit:** Add a protein pair A,123 ↔ B,45 at **6 Å**,
then add another pair A,125 ↔ a selected ligand atom at **8 Å**.
Each card has its own distance field and optional potential. Click **Apply pair constraints** to save both.

**Numbering and supported targets:** Positions start at 1 in each input protein chain,
which can differ from PDB numbering. Protein endpoints select whole residues; specific protein atoms
such as 123/NZ are not supported. Ligand atom names are molecule-specific, so check them when screening
different ligands. Supported limits are **4–20 Å**.

**Strengthen a constraint:** Enable the distance potential to add guidance during sampling.
The potential encourages the limit; it does not guarantee that every prediction satisfies it.
Pocket/contact conditioning and these potentials require **Boltz 2.2 or later**.
""")
        st.markdown("[Boltz2 constraint documentation](https://github.com/jwohlwend/boltz/blob/v2.2.1/docs/prediction.md#constraints)")


def display_binding_pocket_section(protein_sequence: str = None, default_constraints: Any = None) -> List[Dict[str, Any]]:
    """Edit independently configured ligand-to-pocket constraints."""
    applied = normalize_binding_pocket_constraints(
        st.session_state.get("binding_pocket_constraints", default_constraints)
    )
    latest_protein_seq = (
        protein_sequence or st.session_state.get("wt_protein_input")
        or st.session_state.get("protein_input_normal") or ""
    )
    _, _, chains_dict, _ = validate_protein_sequence(latest_protein_seq) if latest_protein_seq else (False, "", {}, "")
    st.subheader("Ligand → binding pocket")
    st.markdown(
        "Guide the **small molecule (chain X)** toward selected protein residues. "
        "The distance limit is measured from the closest ligand atom to the closest atom "
        "in **each listed residue**."
    )
    st.caption(
        "Example: X → A,123 and B,45 at 6 Å asks the ligand to contact each residue within 6 Å. "
        "To constrain A,123 ↔ B,45 directly, use the pair-distance section below."
    )
    st.caption(
        "Add multiple pocket constraints for different residue groups and distance limits. "
        "All applied pocket and pair constraints are used together in the same prediction."
    )
    if "batch_pocket_constraint_rows" not in st.session_state:
        rows = []
        for index, pocket in enumerate(applied):
            rows.append({**pocket, "id": index,
                         "manual_residues": " ".join(f"{c},{r}" for c, r in pocket["contacts"]),
                         "mode": pocket.get("mode", "manual")})
        st.session_state["batch_pocket_constraint_rows"] = rows or [{"id": 0}]
        st.session_state["batch_pocket_constraint_next_id"] = len(rows) or 1
    rows = st.session_state["batch_pocket_constraint_rows"]
    drafts = []
    for index, row in enumerate(rows, 1):
        row_id = row["id"]
        prefix = f"batch_pocket_{row_id}"
        with st.container(border=True):
            st.markdown(f"**Pocket constraint {index} · ligand X → protein residues**")
            selection_col, distance_col = st.columns([2, 1])
            with selection_col:
                mode = st.selectbox(
                    "Choose pocket residues using", ["Residue list", "Sequence search"],
                    index=1 if row.get("mode") == "automated" else 0, key=f"{prefix}_mode",
                )
            with distance_col:
                max_distance = st.number_input(
                    "Ligand-to-pocket limit (Å)", min_value=4.0, max_value=20.0,
                    value=float(row.get("max_distance", 6.0)), step=0.1, key=f"{prefix}_max",
                    help="Closest atom pair between ligand X and each pocket residue; supported range: 4–20 Å.",
                )
            contacts = []
            manual_residues = sequence = ""
            draft_errors = []
            if mode == "Residue list":
                groups_key = f"{prefix}_residue_groups"
                if groups_key not in st.session_state:
                    groups = {}
                    for chain, residue in row.get("contacts", []):
                        groups.setdefault(chain, []).append(str(residue))
                    st.session_state[groups_key] = [
                        {"id": i, "chain": chain, "positions": ", ".join(positions)}
                        for i, (chain, positions) in enumerate(groups.items())
                    ] or [{"id": 0, "chain": "A", "positions": ""}]
                    st.session_state[f"{groups_key}_next_id"] = len(groups) or 1
                groups = st.session_state[groups_key]
                for group in groups:
                    group_prefix = f"{prefix}_group_{group['id']}"
                    chain_col, residue_col, remove_col = st.columns([1, 3, 1], vertical_alignment="bottom")
                    with chain_col:
                        chain_options = sorted(chains_dict) or [chr(i) for i in range(65, 91) if chr(i) != "X"]
                        chain = st.selectbox(
                            "Protein chain", chain_options,
                            index=chain_options.index(group["chain"]) if group["chain"] in chain_options else 0,
                            key=f"{group_prefix}_chain",
                        )
                    with residue_col:
                        positions = st.text_input(
                            "Residue positions", value=group["positions"], placeholder="123, 125, 130-135",
                            key=f"{group_prefix}_positions",
                            help="Enter positions separated by commas or spaces. Ranges such as 130-135 are accepted. Positions start at 1 in the input sequence.",
                        )
                        if positions.strip():
                            manual_residues += f" {chain}:{positions}"
                            try:
                                contacts.extend(parse_pocket_residue_positions(
                                    chain, positions, len(chains_dict[chain]) if chain in chains_dict else None
                                ))
                            except ValueError as exc:
                                draft_errors.append(str(exc))
                    with remove_col:
                        st.button("Remove", key=f"{group_prefix}_remove", disabled=len(groups) == 1, type="tertiary",
                                  on_click=_remove_editor_row, args=(groups_key, group["id"]))
                st.button("Add protein chain", key=f"{prefix}_add_chain", type="tertiary",
                          on_click=_add_pocket_chain, args=(groups_key,))
            else:
                sequence = st.text_input(
                    "Sequence fragments", value=row.get("sequence", ""),
                    placeholder="ACDE FGH", key=f"{prefix}_sequence",
                    help="Separate fragments with spaces or commas. All matched residues share this pocket's distance limit.",
                ).strip().upper()
                fragments = [fragment for fragment in re.split(r"[\s,]+", sequence) if fragment]
                matches = [find_sequence_matches_in_chains(fragment, chains_dict) for fragment in fragments]
                matching_chains = sorted({chain for chains, _ in matches for chain in chains})
                selected_chains = st.multiselect(
                    "Protein chains", matching_chains, default=matching_chains, key=f"{prefix}_chains",
                )
                for fragment, (chains, indices) in zip(fragments, matches):
                    contacts.extend(generate_contacts_from_sequence_matches(fragment, selected_chains, indices, chains_dict))
                    if not chains:
                        draft_errors.append(f"No sequence match for {fragment}. Check the fragment and input protein.")
                if contacts:
                    st.caption("Matched residues: " + ", ".join(dict.fromkeys(f"{c},{r}" for c, r in contacts)))
            for error in draft_errors:
                st.warning(error)
            force = st.toggle(
                "Strengthen this constraint", value=bool(row.get("force", False)),
                key=f"{prefix}_force",
                help="Off: model guidance. On: also apply Boltz's distance potential during sampling. Predicted distances can still vary.",
            )
            if contacts:
                preview_rows = selected_residue_rows(contacts, chains_dict, max_distance)
                st.info(f"Each of the {len(preview_rows)} selected residues should be within {max_distance:g} Å of ligand X.")
                st.markdown("**Confirm selected pocket residues**")
                st.dataframe(preview_rows, hide_index=True, use_container_width=True)
                st.caption("Preview of the current inputs. Click Apply pocket constraints to use this selection.")
            st.button("Remove pocket constraint", key=f"{prefix}_remove", disabled=len(rows) == 1, type="tertiary",
                      on_click=_remove_editor_row, args=("batch_pocket_constraint_rows", row_id))
            drafts.append({"binder": "X", "contacts": contacts, "max_distance": max_distance,
                           "force": force, "mode": "manual" if mode == "Residue list" else "automated",
                           "manual_residues": manual_residues, "sequence": sequence, "_errors": draft_errors})
    add_col, apply_col, clear_col = st.columns(3)
    with add_col:
        if st.button("Add pocket constraint", key="batch_pocket_add"):
            next_id = st.session_state["batch_pocket_constraint_next_id"]
            st.session_state["batch_pocket_constraint_rows"].append({"id": next_id})
            st.session_state["batch_pocket_constraint_next_id"] = next_id + 1
            st.rerun()
    with apply_col:
        apply = st.button("Apply pocket constraints", key="batch_pocket_apply", icon=":material/check_circle:")
    with clear_col:
        if st.button("Clear applied pockets", key="batch_pocket_clear", disabled=not applied, type="tertiary"):
            applied = []
            st.session_state["binding_pocket_constraints"] = []
            st.success("Applied pocket constraints cleared.")
    if apply:
        try:
            populated = [draft for draft in drafts if draft["manual_residues"].strip() or draft["sequence"].strip()]
            errors = [error for draft in populated for error in draft.pop("_errors", [])]
            if errors:
                raise ValueError(" ".join(errors))
            normalized = normalize_binding_pocket_constraints(populated)
            for pocket in normalized:
                for chain, residue in pocket["contacts"]:
                    if not isinstance(residue, int):
                        raise ValueError("Pocket residues must use chain,residue, for example A,123.")
                    if chains_dict and (chain not in chains_dict or residue > len(chains_dict[chain])):
                        raise ValueError(f"Pocket residue {chain},{residue} is outside the input protein sequence.")
        except ValueError as exc:
            st.error(str(exc))
        else:
            applied = normalized
            st.session_state["binding_pocket_constraints"] = applied
            st.success(f"Applied {len(applied)} pocket constraint(s).")
    st.caption("Edits take effect after Apply. Clear removes applied pockets from the next run.")
    display_applied_pocket_constraints(applied)
    return applied


def normalize_distance_constraints(distance_constraints: Any) -> List[Dict[str, Any]]:
    """Validate Boltz residue/atom contact constraints and return YAML-ready values."""
    if distance_constraints is None:
        return []
    if not isinstance(distance_constraints, list):
        raise ValueError("Distance constraints must be a list of constraint objects.")

    def normalize_token(token: Any, label: str) -> List[Any]:
        if not isinstance(token, (list, tuple)) or len(token) != 2:
            raise ValueError(f"{label} must be [CHAIN_ID, RES_IDX or LIGAND_ATOM_NAME].")
        chain_id = str(token[0]).strip()
        if not chain_id:
            raise ValueError(f"{label} needs a chain ID.")
        residue_or_atom = token[1]
        if residue_or_atom is None or isinstance(residue_or_atom, bool):
            raise ValueError(f"{label} residue index must be a positive integer or RES_IDX or LIGAND_ATOM_NAME.")
        if isinstance(residue_or_atom, int):
            if residue_or_atom < 1:
                raise ValueError(f"{label} residue index must be at least 1.")
            normalized_residue_or_atom: Any = residue_or_atom
        else:
            text = str(residue_or_atom).strip()
            if re.fullmatch(r"[1-9]\d*", text):
                normalized_residue_or_atom = int(text)
            elif re.fullmatch(r"[A-Za-z][A-Za-z0-9_'\*]*", text):
                normalized_residue_or_atom = text
            else:
                raise ValueError(
                    f"{label} residue/atom must be a positive residue index (for example 226) "
                    "or a ligand atom name (for example O1). Protein atom selectors such as 77/NZ are not supported by Boltz2."
                )
        return [chain_id, normalized_residue_or_atom]

    normalized = []
    for index, constraint in enumerate(distance_constraints, start=1):
        if not isinstance(constraint, dict):
            raise ValueError(f"Distance constraint {index} must be an object.")
        token1 = normalize_token(constraint.get("token1"), f"Distance constraint {index} endpoint 1")
        token2 = normalize_token(constraint.get("token2"), f"Distance constraint {index} endpoint 2")
        try:
            max_distance = float(constraint.get("max_distance", 6.0))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Distance constraint {index} maximum distance must be a number.") from exc
        if not math.isfinite(max_distance) or not 4.0 <= max_distance <= 20.0:
            raise ValueError(f"Distance constraint {index} maximum distance must be between 4 and 20 Å.")
        force = constraint.get("force", False)
        if not isinstance(force, bool):
            raise ValueError(f"Distance constraint {index} force must be true or false.")
        normalized.append({
            "token1": token1,
            "token2": token2,
            "max_distance": max_distance,
            "force": force,
        })
    seen_pairs = {}
    for index, pair in enumerate(normalized, 1):
        key = frozenset((tuple(pair["token1"]), tuple(pair["token2"])))
        if key in seen_pairs and seen_pairs[key] != pair["max_distance"]:
            raise ValueError(f"Pair constraint {index} repeats the same endpoints with a different limit. Use one limit per pair.")
        seen_pairs[key] = pair["max_distance"]
    return normalized


def display_distance_constraints_section(protein_sequence: str = None, ligand_smiles: Optional[List[Tuple[str, str]]] = None) -> List[Dict[str, Any]]:
    """Edit pair contacts with explicit endpoint types and separate input fields."""
    applied = st.session_state.get("distance_constraints", [])
    if "batch_distance_constraint_rows" not in st.session_state:
        st.session_state["batch_distance_constraint_rows"] = [
            {**constraint, "id": index} for index, constraint in enumerate(applied)
        ] or [{"id": 0}]
        st.session_state["batch_distance_constraint_next_id"] = len(applied) or 1
    latest_sequence = (protein_sequence or st.session_state.get("wt_protein_input")
                       or st.session_state.get("protein_input_normal") or "")
    _, _, protein_chains, _ = validate_protein_sequence(latest_sequence) if latest_sequence else (False, "", {}, "")
    cofactors = st.session_state.get("cofactor_info") or []
    if isinstance(cofactors, dict):
        cofactors = [cofactors]
    ligand_chains = ["X"] + [chr(ord("T") + index) for index, cofactor in enumerate(cofactors[:4])
                             if cofactor.get("smiles") or cofactor.get("ccd")]

    st.subheader("Pair distance · endpoint 1 ↔ endpoint 2")
    st.markdown(
        "Bring **two selected residues or ligand atoms** within a distance limit. "
        "Choose each endpoint below."
    )
    st.caption("Protein residue pairs use the closest atom pair. Two ligand atom endpoints use those exact atoms.")
    rows = st.session_state["batch_distance_constraint_rows"]
    drafts = []
    for index, row in enumerate(rows, 1):
        row_id = row["id"]
        prefix = f"batch_distance_{row_id}"
        with st.container(border=True):
            st.markdown(f"**Pair constraint {index}**")
            endpoint_cols = st.columns(2)
            endpoints = []
            descriptions = []
            for endpoint_index, col in enumerate(endpoint_cols, 1):
                token = row.get(f"token{endpoint_index}", ["A", ""])
                if isinstance(token, str):  # Draft rows from earlier versions.
                    token = token.split(",", 1) if "," in token else ["A", ""]
                chain, target = token
                atom_default = isinstance(target, str) and bool(target) and not target.isdigit()
                endpoint_prefix = f"{prefix}_endpoint{endpoint_index}"
                with col:
                    kind = st.selectbox(
                        f"Endpoint {endpoint_index}", ["Protein residue", "Ligand atom"],
                        index=1 if atom_default else 0, key=f"{endpoint_prefix}_kind",
                    )
                    chain_col, target_col = st.columns([1, 2])
                    with chain_col:
                        options = sorted(protein_chains) or [chr(i) for i in range(65, 91) if chr(i) != "X"]
                        if kind == "Ligand atom":
                            options = ligand_chains
                        chain_value = st.selectbox(
                            "Chain", options, index=options.index(chain) if chain in options else 0,
                            key=f"{endpoint_prefix}_{kind}_chain",
                            format_func=(lambda value: "X · screened ligand" if value == "X" else value),
                        )
                    with target_col:
                        if kind == "Protein residue":
                            residue_default = int(target) if str(target).isdigit() and int(target) > 0 else None
                            target_value = st.number_input(
                                "Residue position", min_value=1, value=residue_default, step=1,
                                placeholder="e.g. 123", key=f"{endpoint_prefix}_residue",
                                help="Position starts at 1 in this chain's input sequence.",
                            )
                        else:
                            preview_ligands = ligand_smiles or []
                            if chain_value != "X":
                                cofactor_index = ord(chain_value) - ord("T")
                                cofactor = cofactors[cofactor_index]
                                preview_ligands = [(f"Co-factor {chain_value}", cofactor["smiles"])] if cofactor.get("smiles") else []
                            if preview_ligands:
                                preview_index = 0
                                if len(preview_ligands) > 1:
                                    preview_index = st.selectbox(
                                        "Show atoms from ligand", range(len(preview_ligands)),
                                        format_func=lambda i, ligands=preview_ligands: ligands[i][0], key=f"{endpoint_prefix}_preview_ligand",
                                    )
                                    st.caption("This atom constraint applies to every screened ligand. Each molecule must contain the selected atom name.")
                                preview_name, preview_smiles = preview_ligands[preview_index]
                                try:
                                    atoms = boltz_ligand_atoms(preview_smiles, affinity=chain_value == "X")
                                except ValueError as exc:
                                    st.error(str(exc))
                                    atoms = ()
                                atom_elements = dict(atoms)
                                atom_options = [None] + list(atom_elements)
                                target_value = st.selectbox(
                                    "Ligand atom", atom_options,
                                    index=atom_options.index(target) if target in atom_options else 0,
                                    format_func=lambda name, elements=atom_elements: "Choose an atom" if name is None else f"{name} ({elements[name]})",
                                    key=f"{endpoint_prefix}_{chain_value}_atom_choice_{preview_index}",
                                    help=f"Names generated using Boltz2's atom naming rules for {preview_name}.",
                                )
                            else:
                                target_value = st.text_input(
                                    "Ligand atom name", value=str(target) if atom_default else "",
                                    placeholder="e.g. O1", key=f"{endpoint_prefix}_atom",
                                    help="Use the exact atom name in Boltz's generated ligand structure or the co-factor CCD definition. Add a ligand SMILES to choose from its atoms.",
                                ).strip()
                    endpoints.append([chain_value, target_value])
                    descriptions.append(
                        f"{'residue' if kind == 'Protein residue' else 'ligand atom'} {chain_value},{target_value}"
                        if target_value else f"endpoint {endpoint_index}"
                    )
            limit_col, potential_col = st.columns(2)
            with limit_col:
                limit = st.number_input(
                    "Distance between these endpoints (Å)", min_value=4.0, max_value=20.0,
                    value=float(row.get("max_distance", 6.0)), step=0.1, key=f"{prefix}_max",
                    help="Maximum distance for this pair only. Supported range: 4–20 Å.",
                )
            with potential_col:
                force = st.toggle(
                    "Strengthen this constraint", value=bool(row.get("force", False)),
                    key=f"{prefix}_force",
                    help="Off: model guidance. On: also apply Boltz's distance potential during sampling. Predicted distances can still vary.",
                )
            st.info(f"{descriptions[0]} ↔ {descriptions[1]} · distance ≤ {limit:g} Å")
            confirmation = []
            for endpoint_number, (chain, target) in enumerate(endpoints, 1):
                if not target:
                    continue
                sequence = protein_chains.get(chain, "")
                residue = sequence[target - 1] if isinstance(target, int) and 1 <= target <= len(sequence) else "Unavailable"
                confirmation.append({"Endpoint": endpoint_number, "Chain": chain,
                                     "Selection": f"Residue {target}" if isinstance(target, int) else f"Ligand atom {target}",
                                     "Residue": residue if isinstance(target, int) else "—", "Pair limit (Å)": limit})
            if confirmation:
                st.markdown("**Confirm selected endpoints**")
                st.dataframe(confirmation, hide_index=True, use_container_width=True)
            st.button("Remove pair constraint", key=f"{prefix}_remove", disabled=len(rows) == 1, type="tertiary",
                      on_click=_remove_editor_row, args=("batch_distance_constraint_rows", row_id))
            drafts.append({"token1": endpoints[0], "token2": endpoints[1], "max_distance": limit, "force": force})
    add_col, apply_col, clear_col = st.columns(3)
    with add_col:
        if st.button("Add pair constraint", key="batch_distance_add"):
            next_id = st.session_state["batch_distance_constraint_next_id"]
            st.session_state["batch_distance_constraint_rows"].append({"id": next_id})
            st.session_state["batch_distance_constraint_next_id"] = next_id + 1
            st.rerun()
    with apply_col:
        apply = st.button("Apply pair constraints", key="batch_distance_apply", icon=":material/check_circle:")
    with clear_col:
        if st.button("Clear applied pairs", key="batch_distance_clear", disabled=not applied, type="tertiary"):
            applied = []
            st.session_state["distance_constraints"] = []
            st.success("Applied pair constraints cleared.")
    if apply:
        try:
            populated = []
            for index, draft in enumerate(drafts, 1):
                first, second = draft["token1"], draft["token2"]
                if not first[1] and not second[1]:
                    continue
                if not first[1] or not second[1]:
                    raise ValueError(f"Pair constraint {index}: choose both endpoints before applying.")
                for chain, target in (first, second):
                    if chain in ligand_chains and isinstance(target, str) and target.isdigit():
                        raise ValueError(f"Ligand chain {chain}: enter an atom name such as O1, not a residue number.")
                    if chain in protein_chains and isinstance(target, int) and target > len(protein_chains[chain]):
                        raise ValueError(f"Residue {chain},{target} is outside the input protein sequence.")
                if first == second:
                    raise ValueError(f"Pair constraint {index}: select two different endpoints.")
                populated.append(draft)
            normalized = normalize_distance_constraints(populated)
        except ValueError as exc:
            st.error(str(exc))
        else:
            applied = normalized
            st.session_state["distance_constraints"] = applied
            st.success(f"Applied {len(applied)} pair constraint(s).")
    st.caption("Edits take effect after Apply. All applied pockets and pairs are used together. Clear removes applied pairs from the next run.")
    if applied:
        st.markdown("**Applied pair-distance constraints**")
        for constraint in applied:
            first, second = constraint["token1"], constraint["token2"]
            enforcement = "potential enabled" if constraint.get("force") else "model guidance"
            st.write(f"{first[0]},{first[1]} ↔ {second[0]},{second[1]} ≤ {constraint['max_distance']:g} Å · {enforcement}")
    return applied


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
        mutation_parts = [s.strip() for s in re.split(r"[-/+]", mutant_str)]
        
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
        
        # Preserve one entry per comma-delimited mutant so invalid entries do
        # not shift subsequent mutations onto the wrong label.
        mutants.append(mutations)
    
    return mutants


def suggest_chain_start_from_mutations(mutation_text: str, chain_sequence: str) -> Dict[str, Any]:
    """Infer possible 1-based chain starts from mutation WT residues.

    A start is returned only when one offset makes every entered mutation's
    wild-type residue agree with the supplied chain sequence. Multiple exact
    offsets are reported as ambiguous so callers can ask the user to choose.
    """
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    entries: List[Tuple[str, int]] = []
    invalid_entries: List[str] = []

    for mutant in re.split(r"[,;]", str(mutation_text or "")):
        mutant = mutant.strip()
        if not mutant:
            continue
        for part in re.split(r"[/+\-]", mutant):
            part = part.strip().upper()
            if not part:
                invalid_entries.append(part)
                continue
            match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", part)
            if not match:
                invalid_entries.append(part)
                continue
            wt_residue, position_text, _new_residue = match.groups()
            if wt_residue not in amino_acids or _new_residue not in amino_acids:
                invalid_entries.append(part)
                continue
            try:
                position = int(position_text)
            except (TypeError, ValueError, OverflowError):
                invalid_entries.append(part)
                continue
            if position < 1:
                invalid_entries.append(part)
                continue
            entries.append((wt_residue, position))

    # Repeated alternatives at the same WT position provide no extra evidence.
    entries = list(dict.fromkeys(entries))
    if invalid_entries:
        return {
            "status": "invalid_mutations",
            "candidates": [],
            "mutation_count": len(entries),
            "invalid_entries": invalid_entries,
        }
    if not entries:
        return {
            "status": "no_mutations",
            "candidates": [],
            "mutation_count": 0,
            "invalid_entries": [],
        }

    sequence = re.sub(r"\s+", "", str(chain_sequence or "")).upper()
    if not sequence or not re.fullmatch(r"[A-Z]+", sequence):
        return {
            "status": "invalid_sequence",
            "candidates": [],
            "mutation_count": len(entries),
            "invalid_entries": [],
        }

    possible_starts: Optional[set] = None
    for wt_residue, position in entries:
        mutation_starts = {
            position - index
            for index, residue in enumerate(sequence)
            if residue == wt_residue and position - index >= 1
        }
        possible_starts = (
            mutation_starts
            if possible_starts is None
            else possible_starts.intersection(mutation_starts)
        )
        if not possible_starts:
            break

    candidates = sorted(possible_starts or [])
    return {
        "status": "unique" if len(candidates) == 1 else "ambiguous" if candidates else "no_match",
        "candidates": candidates,
        "mutation_count": len(entries),
        "invalid_entries": [],
    }

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
        mutation_parts = [s.strip() for s in re.split(r"[-/+]", mutant_str)]
        
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
        mutation_parts = [s.strip() for s in re.split(r"[-/+]", mutation_str)]
        
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
            mapping[i + 1] = i + 1
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
            mapping[canonical_pos] = i + 1

    return mapping


def map_canonical_to_input_position(canonical_pos: int, mapping: Dict[int, int], chain_start: int) -> Optional[int]:
    """Map canonical position to input sequence position."""
    if canonical_pos in mapping:
        return mapping[canonical_pos]

    # Fallback: assume 1:1 local sequence mapping if not found
    return canonical_pos if canonical_pos > 0 else None


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
                                                    'position_in_user_sequence': input_pos,
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
                                            input_mutation_pos = chain_start + input_pos - 1
                                            # Categorize the change
                                            change_type = categorize_amino_acid_change(wt_residue, mut_residue)
                                            severity = predict_functional_severity(wt_residue, mut_residue, change_type)
                                            structural_impact = predict_structural_impact(change_type, severity)

                                            mutation_data = {
                                                'mutation_code': f"{wt_residue}{input_mutation_pos}{mut_residue}",
                                                'canonical_code': f"{wt_residue}{canonical_pos}{mut_residue}",
                                                'input_position': input_mutation_pos,
                                                'position_in_user_sequence': input_pos,
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
        with st.expander("AlphaMissense Prioritization (Optional)", expanded=False):
            st.caption(
                "AlphaMissense controls appear here once a wild-type sequence is provided "
                "and mutation discovery is available."
            )
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

    # Optional AlphaMissense integration controls
    am_lookup = get_default_alphamissense_lookup() if get_default_alphamissense_lookup is not None else None
    am_backend_available = bool(am_lookup is not None and getattr(am_lookup, "available", False))
    am_backend_reason = (
        getattr(am_lookup, "status_reason", "unavailable")
        if am_lookup is not None
        else "module_import_failed"
    )

    with st.expander("AlphaMissense Prioritization (Optional)", expanded=False):
        if am_backend_available:
            st.caption("Local AlphaMissense index detected. Variant-level pathogenicity scores are available.")
        else:
            st.caption(
                "AlphaMissense lookup is unavailable. Install `pysam` and keep the `.tsv.gz` + `.tbi` files "
                "in `application/alphamissense/` to enable indexed lookups."
            )
            st.caption(f"Status: `{am_backend_reason}`")

        use_alphamissense_prior = st.toggle(
            "Use AlphaMissense-Informed Prioritization",
            value=st.session_state.get("use_alphamissense_prior", True),
            key="use_alphamissense_prior",
            disabled=not am_backend_available,
            help=(
                "Ranks discovered variants by a combined score: AlphaMissense pathogenicity + existing "
                "severity/evidence annotations."
            ),
        )
        use_priority_sort = st.toggle(
            "Sort by Priority Score",
            value=st.session_state.get("use_alphamissense_priority_sort", True),
            key="use_alphamissense_priority_sort",
            disabled=not use_alphamissense_prior,
            help="If OFF, mutations are shown in numeric residue order.",
        )
        apply_am_filters = st.toggle(
            "Filter by AlphaMissense Criteria",
            value=st.session_state.get("use_alphamissense_filters", False),
            key="use_alphamissense_filters",
            disabled=not use_alphamissense_prior,
            help="Optionally keep only mutations above a score threshold and selected classes.",
        )
        am_min_score = st.slider(
            "Minimum AlphaMissense pathogenicity score",
            min_value=0.0,
            max_value=1.0,
            value=float(st.session_state.get("alphamissense_min_score", 0.60)),
            step=0.01,
            key="alphamissense_min_score",
            disabled=not (use_alphamissense_prior and apply_am_filters),
        )
        am_allowed_classes = st.multiselect(
            "Allowed AlphaMissense classes",
            options=["pathogenic", "ambiguous", "benign"],
            default=st.session_state.get("alphamissense_allowed_classes", ["pathogenic", "ambiguous"]),
            key="alphamissense_allowed_classes",
            disabled=not (use_alphamissense_prior and apply_am_filters),
            help="Class labels as provided by AlphaMissense.",
        )
        prioritized_copy_count = st.slider(
            "Top prioritized mutations for quick copy",
            min_value=5,
            max_value=100,
            value=int(st.session_state.get("alphamissense_copy_top_n", 20)),
            step=1,
            key="alphamissense_copy_top_n",
            disabled=not use_alphamissense_prior,
        )

    # Enforce backend availability even if stale session state preserved an enabled toggle.
    use_alphamissense_prior = bool(use_alphamissense_prior and am_backend_available)
    if not use_alphamissense_prior:
        use_priority_sort = False
        apply_am_filters = False

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
    if 'discovered_mutations_raw' not in st.session_state:
        st.session_state.discovered_mutations_raw = None
    if 'discovery_sequence_hash' not in st.session_state:
        st.session_state.discovery_sequence_hash = None

    # Create a hash of the current sequence to detect changes
    import hashlib
    current_seq_hash = hashlib.md5(first_chain_seq.encode()).hexdigest()

    # Clear stored mutations if sequence changed
    if st.session_state.discovery_sequence_hash != current_seq_hash:
        st.session_state.discovered_mutations = None
        st.session_state.discovered_mutations_raw = None
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
            st.session_state.discovered_mutations_raw = detailed_mutations
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
            st.session_state.discovered_mutations_raw = detailed_mutations
            st.session_state.discovered_mutations = detailed_mutations
            st.session_state.discovery_sequence_hash = current_seq_hash

    # Display mutations if they exist in session state (persists across reruns)
    if st.session_state.discovered_mutations_raw:
        detailed_mutations = [dict(item) for item in st.session_state.discovered_mutations_raw]
        am_summary: Dict[str, Any] = {
            "available": False,
            "reason": "disabled",
            "matched": 0,
            "total": len(detailed_mutations),
            "returned": len(detailed_mutations),
        }
        if use_alphamissense_prior and annotate_mutations_with_alphamissense is not None:
            am_min_filter = float(am_min_score) if apply_am_filters else None
            am_classes_filter = am_allowed_classes if apply_am_filters else None
            detailed_mutations, am_summary = annotate_mutations_with_alphamissense(
                detailed_mutations,
                use_prioritization=use_priority_sort,
                min_am_score=am_min_filter,
                allowed_am_classes=am_classes_filter,
                top_n=None,
            )
            st.session_state["mutation_discovery_priority_map"] = {
                str(item.get("mutation_code", "")): float(item.get("priority_score", 0.0))
                for item in detailed_mutations
                if item.get("mutation_code")
            }
            st.session_state["mutation_discovery_alpha_summary"] = am_summary
            if am_summary.get("available"):
                st.caption(
                    f"AlphaMissense matches: {am_summary.get('matched', 0)}/{am_summary.get('total', 0)}"
                    f" variants (returned {am_summary.get('returned', 0)} after filters)."
                )
            else:
                st.caption(f"AlphaMissense unavailable: {am_summary.get('reason')}")
        else:
            st.session_state["mutation_discovery_priority_map"] = {}
            st.session_state["mutation_discovery_alpha_summary"] = am_summary

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
        if use_alphamissense_prior and use_priority_sort:
            filtered_mutations.sort(
                key=lambda x: (
                    float(x.get("priority_score", 0.0)),
                    float(x.get("evidence_count", 0.0)),
                ),
                reverse=True,
            )
        else:
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
                "AlphaMissense Score": (
                    f"{float(m.get('am_pathogenicity')):.3f}"
                    if m.get('am_pathogenicity') is not None
                    else '-'
                ),
                "AlphaMissense Class": m.get('am_class', '-') if m.get('am_class') else '-',
                "Priority Score": f"{float(m.get('priority_score', 0.0)):.3f}",
                "Priority Tier": m.get('priority_tier', 'LOW'),
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
                "AlphaMissense Score": st.column_config.TextColumn(
                    "AlphaMissense Score",
                    width="small",
                    help="AlphaMissense pathogenicity score (0 to 1). Higher is more likely pathogenic.",
                ),
                "AlphaMissense Class": st.column_config.TextColumn(
                    "AlphaMissense Class",
                    width="small",
                    help="AlphaMissense class label: pathogenic / ambiguous / benign.",
                ),
                "Priority Score": st.column_config.TextColumn(
                    "Priority Score",
                    width="small",
                    help="Combined prioritization score using AlphaMissense, SIFT, PolyPhen, clinical significance, and legacy severity/evidence.",
                ),
                "Priority Tier": st.column_config.TextColumn(
                    "Priority",
                    width="small",
                    help="Derived triage tier from the priority score.",
                ),
                "Clinical Significance": st.column_config.TextColumn("Clinical Significance", width="medium"),
                "References": st.column_config.TextColumn("References", width="medium")
            }
        )

        # Format mutations for input (display below table)
        st.markdown("##### :material/content_copy: Copy-Pastable Mutation List")
        formatted_mutations = format_mutations_for_input(filtered_mutations)
        st.code(formatted_mutations, language=None)
        if use_alphamissense_prior:
            prioritized_for_copy = sorted(
                filtered_mutations,
                key=lambda x: (
                    float(x.get("priority_score", 0.0)),
                    float(x.get("evidence_count", 0.0)),
                ),
                reverse=True,
            )[: max(1, int(prioritized_copy_count))]
            st.markdown("##### :material/target: Priority-Selected Mutation List")
            st.code(format_mutations_for_input(prioritized_for_copy), language=None)

        # Column explanations
        st.markdown("##### :material/info: Column Descriptions")
        st.markdown("""
        - **Mutation (Your Seq)**: Mutation using your sequence numbering (with applied offset). **Use this for mutations input.**
        - **Canonical**: Mutation in standard UniProt numbering - for reference lookup
        - **SIFT**: Predicts tolerance - deleterious if score ≤0.05 (may be empty if not available in database)
        - **PolyPhen**: Predicts damage - probably/possibly damaging or benign (may be empty if not available)
        - **AlphaMissense Score / AlphaMissense Class**: AlphaMissense pathogenicity evidence from local indexed dataset
        - **Priority Score / Priority**: Combined triage score incorporating AlphaMissense, SIFT, PolyPhen, clinical significance, and legacy severity
        - **Clinical Significance**: Known clinical impact from databases
        - **References**: References (PMID or database IDs) supporting this variant

        **Numbering Adjustment:** Use the "Adjust Numbering & Filters" controls above to set a position offset if your sequence
        numbering differs from the canonical UniProt positions. The suggested offset aligns your sequence with the full protein.
        Apply filters to show only mutations in a specific residue range.
        """)

        return formatted_mutations

    return None
