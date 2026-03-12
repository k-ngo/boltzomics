#!/usr/bin/env python3
"""
Mutation Analysis Module for BoltzOmics

This module provides analysis of mutation effects on protein structure and
drug binding, including:

1. Binding site proximity analysis - Is the mutation near the drug?
2. Mutation property analysis - Physicochemical changes from mutation
3. Residue numbering utilities - Map between canonical/input/PDB numbering

Author: BoltzOmics Team
"""

import re
import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any
from enum import Enum
import numpy as np

logger = logging.getLogger(__name__)


# ============================================================================
# Amino Acid Properties
# ============================================================================

# Standard amino acid properties
AA_PROPERTIES = {
    # Amino acid: (charge at pH 7, hydrophobicity, molecular weight, volume Å³, is_aromatic)
    'A': (0, 1.8, 89, 88.6, False),    # Alanine
    'R': (1, -4.5, 174, 173.4, False),  # Arginine
    'N': (0, -3.5, 132, 114.1, False),  # Asparagine
    'D': (-1, -3.5, 133, 111.1, False), # Aspartic acid
    'C': (0, 2.5, 121, 108.5, False),   # Cysteine
    'E': (-1, -3.5, 147, 138.4, False), # Glutamic acid
    'Q': (0, -3.5, 146, 143.8, False),  # Glutamine
    'G': (0, -0.4, 75, 60.1, False),    # Glycine
    'H': (0.1, -3.2, 155, 153.2, True), # Histidine
    'I': (0, 4.5, 131, 166.7, False),   # Isoleucine
    'L': (0, 3.8, 131, 166.7, False),   # Leucine
    'K': (1, -3.9, 146, 168.6, False),  # Lysine
    'M': (0, 1.9, 149, 162.9, False),   # Methionine
    'F': (0, 2.8, 165, 189.9, True),    # Phenylalanine
    'P': (0, -1.6, 115, 112.7, False),  # Proline
    'S': (0, -0.8, 105, 89.0, False),   # Serine
    'T': (0, -0.7, 119, 116.1, False),  # Threonine
    'W': (0, -0.9, 204, 227.8, True),   # Tryptophan
    'Y': (0, -1.3, 181, 193.6, True),   # Tyrosine
    'V': (0, 4.2, 117, 140.0, False),   # Valine
}

# Amino acid categories
AA_CATEGORIES = {
    'hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'],
    'polar': ['S', 'T', 'N', 'Q', 'C', 'Y'],
    'positive': ['K', 'R', 'H'],
    'negative': ['D', 'E'],
    'aromatic': ['F', 'Y', 'W', 'H'],
    'small': ['G', 'A', 'S'],
    'large': ['F', 'Y', 'W', 'R', 'K'],
}


# ============================================================================
# Data Classes
# ============================================================================

class ProximityCategory(Enum):
    """Categories for mutation proximity to binding site."""
    DIRECT_CONTACT = "direct_contact"      # < 4 Å
    FIRST_SHELL = "first_shell"            # 4-6 Å
    SECOND_SHELL = "second_shell"          # 6-10 Å
    DISTANT = "distant"                    # > 10 Å


@dataclass
class MutationPosition:
    """Represents a mutation with different numbering schemes."""
    canonical_position: int      # UniProt/canonical numbering (e.g., 324)
    input_position: int          # User's input sequence numbering
    pdb_position: int            # Position in PDB file (1-indexed from sequence start)
    wt_residue: str              # Wild-type amino acid
    mut_residue: str             # Mutant amino acid
    chain_start: int             # Starting residue number of user's sequence

    @property
    def mutation_code(self) -> str:
        """Return mutation in standard format (e.g., G324A)."""
        return f"{self.wt_residue}{self.canonical_position}{self.mut_residue}"

    @property
    def pdb_mutation_code(self) -> str:
        """Return mutation using PDB numbering."""
        return f"{self.wt_residue}{self.pdb_position}{self.mut_residue}"


@dataclass
class ProximityAnalysis:
    """Results of binding site proximity analysis."""
    mutation: MutationPosition
    min_distance_to_ligand: float  # Angstroms
    category: ProximityCategory
    nearby_ligand_atoms: List[str]
    contact_residues_nearby: List[str]  # Other residues near mutation
    is_in_binding_pocket: bool
    interpretation: str


@dataclass
class MutationPropertyChange:
    """Analysis of physicochemical property changes from mutation."""
    mutation: MutationPosition
    charge_change: float           # e.g., +1 for K→R is 0, K→E is -2
    hydrophobicity_change: float   # Kyte-Doolittle scale change
    volume_change: float           # Å³
    weight_change: float           # Daltons
    aromatic_change: str           # "gained", "lost", "unchanged"
    category_changes: List[str]    # e.g., ["hydrophobic→polar", "small→large"]
    severity_score: float          # 0-1, higher = more disruptive
    interpretation: str


# ============================================================================
# Residue Numbering Utilities
# ============================================================================

def parse_mutation(mutation_str: str) -> Tuple[str, int, str]:
    """
    Parse mutation string into components.

    Args:
        mutation_str: Mutation like "G324A" or "Gly324Ala"

    Returns:
        Tuple of (wt_residue, position, mut_residue)
    """
    # Try single-letter format first: G324A
    match = re.match(r'^([A-Z])(\d+)([A-Z])$', mutation_str.upper())
    if match:
        return match.group(1), int(match.group(2)), match.group(3)

    # Try three-letter format: Gly324Ala
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    match = re.match(r'^([A-Za-z]{3})(\d+)([A-Za-z]{3})$', mutation_str)
    if match:
        wt = three_to_one.get(match.group(1).upper())
        mut = three_to_one.get(match.group(3).upper())
        if wt and mut:
            return wt, int(match.group(2)), mut

    raise ValueError(f"Could not parse mutation: {mutation_str}")


def canonical_to_pdb_position(
    canonical_pos: int,
    chain_start: int,
    sequence_length: int
) -> Optional[int]:
    """
    Convert canonical (UniProt) position to PDB residue number.

    Args:
        canonical_pos: Position in canonical numbering
        chain_start: Starting residue number of user's sequence (canonical)
        sequence_length: Length of the input sequence

    Returns:
        PDB residue number (1-indexed) or None if out of range
    """
    # PDB position = canonical_pos - chain_start + 1
    pdb_pos = canonical_pos - chain_start + 1

    if 1 <= pdb_pos <= sequence_length:
        return pdb_pos
    return None


def pdb_to_canonical_position(
    pdb_pos: int,
    chain_start: int
) -> int:
    """
    Convert PDB residue number to canonical (UniProt) position.

    Args:
        pdb_pos: Position in PDB file (1-indexed)
        chain_start: Starting residue number of user's sequence (canonical)

    Returns:
        Canonical position
    """
    return pdb_pos + chain_start - 1


def create_mutation_position(
    mutation_str: str,
    chain_start: int,
    sequence_length: int,
    is_canonical: bool = True
) -> MutationPosition:
    """
    Create a MutationPosition object with all numbering schemes.

    Args:
        mutation_str: Mutation like "G324A"
        chain_start: Starting residue of user's sequence in canonical numbering
        sequence_length: Length of the input sequence
        is_canonical: If True, mutation uses canonical numbering; if False, uses input numbering

    Returns:
        MutationPosition with all numbering schemes
    """
    wt, pos, mut = parse_mutation(mutation_str)

    if is_canonical:
        canonical_pos = pos
        pdb_pos = canonical_to_pdb_position(pos, chain_start, sequence_length)
        input_pos = pos  # Assuming input uses canonical
    else:
        # Mutation is in input/PDB numbering
        pdb_pos = pos
        canonical_pos = pdb_to_canonical_position(pos, chain_start)
        input_pos = pos

    return MutationPosition(
        canonical_position=canonical_pos,
        input_position=input_pos,
        pdb_position=pdb_pos if pdb_pos else pos,
        wt_residue=wt,
        mut_residue=mut,
        chain_start=chain_start
    )


# ============================================================================
# Binding Site Proximity Analysis
# ============================================================================

def analyze_binding_proximity(
    pdb_path: str,
    mutation_pdb_position: int,
    protein_chain: str = "A",
    ligand_chain: str = None,
    contact_threshold: float = 4.0,
    first_shell: float = 6.0,
    second_shell: float = 10.0
) -> ProximityAnalysis:
    """
    Analyze how close a mutation is to the ligand binding site.

    Args:
        pdb_path: Path to PDB file
        mutation_pdb_position: Residue number in PDB file
        protein_chain: Chain ID for protein
        ligand_chain: Chain ID for ligand (auto-detected if None)
        contact_threshold: Distance for direct contact (Å)
        first_shell: Distance for first coordination shell (Å)
        second_shell: Distance for second coordination shell (Å)

    Returns:
        ProximityAnalysis with distance and interpretation
    """
    # Auto-detect chains if needed
    if ligand_chain is None:
        protein_chain, ligand_chain = _detect_chains(pdb_path)

    # Parse PDB
    mutation_atoms, ligand_atoms, other_residues = _parse_pdb_for_proximity(
        pdb_path, mutation_pdb_position, protein_chain, ligand_chain
    )

    if not mutation_atoms:
        return ProximityAnalysis(
            mutation=None,
            min_distance_to_ligand=float('inf'),
            category=ProximityCategory.DISTANT,
            nearby_ligand_atoms=[],
            contact_residues_nearby=[],
            is_in_binding_pocket=False,
            interpretation=f"Mutation at position {mutation_pdb_position} not found in structure"
        )

    if not ligand_atoms:
        return ProximityAnalysis(
            mutation=None,
            min_distance_to_ligand=float('inf'),
            category=ProximityCategory.DISTANT,
            nearby_ligand_atoms=[],
            contact_residues_nearby=[],
            is_in_binding_pocket=False,
            interpretation="No ligand found in structure"
        )

    # Calculate minimum distance to ligand
    min_distance = float('inf')
    nearby_atoms = []

    for m_atom in mutation_atoms:
        for l_atom in ligand_atoms:
            dist = np.linalg.norm(
                np.array(m_atom['coord']) - np.array(l_atom['coord'])
            )
            if dist < min_distance:
                min_distance = dist
            if dist < second_shell:
                nearby_atoms.append(f"{l_atom['name']}:{dist:.1f}Å")

    # Determine category
    if min_distance < contact_threshold:
        category = ProximityCategory.DIRECT_CONTACT
        interpretation = (
            f"CRITICAL: Mutation is in direct contact with ligand ({min_distance:.1f}Å). "
            f"High likelihood of affecting binding affinity."
        )
    elif min_distance < first_shell:
        category = ProximityCategory.FIRST_SHELL
        interpretation = (
            f"IMPORTANT: Mutation is in first coordination shell ({min_distance:.1f}Å). "
            f"May affect binding through direct interactions or local conformational changes."
        )
    elif min_distance < second_shell:
        category = ProximityCategory.SECOND_SHELL
        interpretation = (
            f"MODERATE: Mutation is in second shell ({min_distance:.1f}Å). "
            f"Could indirectly affect binding through conformational propagation."
        )
    else:
        category = ProximityCategory.DISTANT
        interpretation = (
            f"LOW IMPACT EXPECTED: Mutation is distant from binding site ({min_distance:.1f}Å). "
            f"Unlikely to directly affect ligand binding."
        )

    # Check if mutation is in binding pocket (contacting residues)
    is_in_pocket = min_distance < first_shell

    # Find nearby residues
    nearby_residues = _find_nearby_residues(
        mutation_atoms, other_residues, cutoff=6.0
    )

    return ProximityAnalysis(
        mutation=None,  # Caller should set this
        min_distance_to_ligand=min_distance,
        category=category,
        nearby_ligand_atoms=nearby_atoms[:5],  # Top 5
        contact_residues_nearby=nearby_residues[:10],
        is_in_binding_pocket=is_in_pocket,
        interpretation=interpretation
    )


def _detect_chains(pdb_path: str) -> Tuple[str, str]:
    """Auto-detect protein and ligand chain IDs."""
    chain_atom_counts = {}
    chain_hetatm_counts = {}

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                chain_atom_counts[chain] = chain_atom_counts.get(chain, 0) + 1
            elif line.startswith('HETATM'):
                chain = line[21]
                resname = line[17:20].strip()
                if resname not in ['HOH', 'WAT', 'NA', 'CL', 'MG', 'ZN', 'CA']:
                    chain_hetatm_counts[chain] = chain_hetatm_counts.get(chain, 0) + 1

    protein_chain = max(chain_atom_counts, key=chain_atom_counts.get) if chain_atom_counts else 'A'

    if chain_hetatm_counts:
        non_protein = {k: v for k, v in chain_hetatm_counts.items() if k != protein_chain}
        if non_protein:
            ligand_chain = max(non_protein, key=non_protein.get)
        else:
            ligand_chain = max(chain_hetatm_counts, key=chain_hetatm_counts.get)
    else:
        other = {k: v for k, v in chain_atom_counts.items() if k != protein_chain}
        ligand_chain = max(other, key=other.get) if other else 'B'

    return protein_chain, ligand_chain


def _parse_pdb_for_proximity(
    pdb_path: str,
    mutation_resid: int,
    protein_chain: str,
    ligand_chain: str
) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    """Parse PDB and extract mutation, ligand, and other residue atoms."""
    mutation_atoms = []
    ligand_atoms = []
    other_residues = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'name': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': line[21],
                    'resid': int(line[22:26].strip()),
                    'coord': [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ]
                }

                if atom['chain'] == protein_chain:
                    if atom['resid'] == mutation_resid:
                        mutation_atoms.append(atom)
                    else:
                        other_residues.append(atom)
                elif atom['chain'] == ligand_chain:
                    ligand_atoms.append(atom)

    return mutation_atoms, ligand_atoms, other_residues


def _find_nearby_residues(
    mutation_atoms: List[Dict],
    other_atoms: List[Dict],
    cutoff: float
) -> List[str]:
    """Find residues near the mutation site."""
    nearby = set()

    for m_atom in mutation_atoms:
        for o_atom in other_atoms:
            dist = np.linalg.norm(
                np.array(m_atom['coord']) - np.array(o_atom['coord'])
            )
            if dist < cutoff:
                nearby.add(f"{o_atom['resname']}{o_atom['resid']}")

    return sorted(list(nearby))


# ============================================================================
# Mutation Property Analysis
# ============================================================================

def analyze_mutation_properties(
    wt_residue: str,
    mut_residue: str,
    position: int = 0
) -> MutationPropertyChange:
    """
    Analyze the physicochemical property changes from a mutation.

    Args:
        wt_residue: Wild-type amino acid (single letter)
        mut_residue: Mutant amino acid (single letter)
        position: Residue position (for reporting)

    Returns:
        MutationPropertyChange with detailed property analysis
    """
    wt = wt_residue.upper()
    mut = mut_residue.upper()

    if wt not in AA_PROPERTIES or mut not in AA_PROPERTIES:
        return MutationPropertyChange(
            mutation=None,
            charge_change=0,
            hydrophobicity_change=0,
            volume_change=0,
            weight_change=0,
            aromatic_change="unknown",
            category_changes=[],
            severity_score=0.5,
            interpretation=f"Unknown amino acid: {wt} or {mut}"
        )

    wt_props = AA_PROPERTIES[wt]
    mut_props = AA_PROPERTIES[mut]

    # Calculate changes
    charge_change = mut_props[0] - wt_props[0]
    hydro_change = mut_props[1] - wt_props[1]
    volume_change = mut_props[3] - wt_props[3]
    weight_change = mut_props[2] - wt_props[2]

    # Aromatic change
    if wt_props[4] and not mut_props[4]:
        aromatic_change = "lost"
    elif not wt_props[4] and mut_props[4]:
        aromatic_change = "gained"
    else:
        aromatic_change = "unchanged"

    # Category changes
    category_changes = []
    for cat_name, cat_residues in AA_CATEGORIES.items():
        wt_in = wt in cat_residues
        mut_in = mut in cat_residues
        if wt_in and not mut_in:
            category_changes.append(f"lost_{cat_name}")
        elif not wt_in and mut_in:
            category_changes.append(f"gained_{cat_name}")

    # Calculate severity score (0-1)
    severity = 0.0

    # Charge changes are severe
    if abs(charge_change) >= 2:
        severity += 0.4
    elif abs(charge_change) >= 1:
        severity += 0.25

    # Large hydrophobicity changes
    if abs(hydro_change) >= 4:
        severity += 0.25
    elif abs(hydro_change) >= 2:
        severity += 0.15

    # Volume changes
    if abs(volume_change) >= 80:
        severity += 0.2
    elif abs(volume_change) >= 40:
        severity += 0.1

    # Aromatic changes
    if aromatic_change in ["lost", "gained"]:
        severity += 0.15

    severity = min(1.0, severity)

    # Generate interpretation
    interpretation_parts = []

    if charge_change > 0:
        interpretation_parts.append(f"introduces positive charge (+{charge_change})")
    elif charge_change < 0:
        interpretation_parts.append(f"introduces negative charge ({charge_change})")

    if hydro_change > 2:
        interpretation_parts.append("increases hydrophobicity")
    elif hydro_change < -2:
        interpretation_parts.append("decreases hydrophobicity")

    if volume_change > 40:
        interpretation_parts.append(f"larger side chain (+{volume_change:.0f}Å³)")
    elif volume_change < -40:
        interpretation_parts.append(f"smaller side chain ({volume_change:.0f}Å³)")

    if aromatic_change == "lost":
        interpretation_parts.append("loses aromatic ring (may affect π-stacking)")
    elif aromatic_change == "gained":
        interpretation_parts.append("gains aromatic ring")

    if severity >= 0.5:
        severity_desc = "HIGH IMPACT"
    elif severity >= 0.25:
        severity_desc = "MODERATE IMPACT"
    else:
        severity_desc = "LOW IMPACT"

    if interpretation_parts:
        interpretation = f"{severity_desc}: {wt}{position}{mut} - " + ", ".join(interpretation_parts)
    else:
        interpretation = f"{severity_desc}: {wt}{position}{mut} - Conservative substitution"

    return MutationPropertyChange(
        mutation=None,
        charge_change=charge_change,
        hydrophobicity_change=hydro_change,
        volume_change=volume_change,
        weight_change=weight_change,
        aromatic_change=aromatic_change,
        category_changes=category_changes,
        severity_score=severity,
        interpretation=interpretation
    )


# ============================================================================
# Combined Analysis
# ============================================================================

def analyze_mutation_impact(
    pdb_path: str,
    mutation_str: str,
    chain_start: int = 1,
    sequence_length: int = None,
    is_canonical: bool = True
) -> Dict[str, Any]:
    """
    Comprehensive mutation impact analysis combining proximity and properties.

    Args:
        pdb_path: Path to predicted structure PDB
        mutation_str: Mutation like "G324A"
        chain_start: Starting residue of input sequence in canonical numbering
        sequence_length: Length of input sequence (for validation)
        is_canonical: Whether mutation uses canonical numbering

    Returns:
        Dictionary with comprehensive analysis
    """
    # Parse mutation
    wt, pos, mut = parse_mutation(mutation_str)

    # Calculate PDB position
    if is_canonical:
        pdb_position = pos - chain_start + 1
    else:
        pdb_position = pos

    # Validate position if sequence length provided
    if sequence_length and (pdb_position < 1 or pdb_position > sequence_length):
        return {
            "error": f"Mutation position {pos} is outside sequence range",
            "mutation": mutation_str,
            "pdb_position": pdb_position,
            "valid_range": f"1-{sequence_length}"
        }

    # Property analysis
    property_analysis = analyze_mutation_properties(wt, mut, pos)

    # Proximity analysis
    proximity_analysis = analyze_binding_proximity(pdb_path, pdb_position)

    # Combined risk assessment
    proximity_risk = {
        ProximityCategory.DIRECT_CONTACT: 1.0,
        ProximityCategory.FIRST_SHELL: 0.7,
        ProximityCategory.SECOND_SHELL: 0.4,
        ProximityCategory.DISTANT: 0.1
    }

    combined_risk = (
        0.5 * property_analysis.severity_score +
        0.5 * proximity_risk.get(proximity_analysis.category, 0.1)
    )

    # Risk interpretation
    if combined_risk >= 0.7:
        risk_level = "HIGH"
        risk_interpretation = "Mutation likely to significantly affect drug binding"
    elif combined_risk >= 0.4:
        risk_level = "MODERATE"
        risk_interpretation = "Mutation may affect drug binding"
    else:
        risk_level = "LOW"
        risk_interpretation = "Mutation unlikely to significantly affect drug binding"

    return {
        "mutation": mutation_str,
        "canonical_position": pos if is_canonical else pdb_position + chain_start - 1,
        "pdb_position": pdb_position,
        "chain_start": chain_start,

        "proximity": {
            "distance_to_ligand": proximity_analysis.min_distance_to_ligand,
            "category": proximity_analysis.category.value,
            "is_in_binding_pocket": proximity_analysis.is_in_binding_pocket,
            "nearby_ligand_atoms": proximity_analysis.nearby_ligand_atoms,
            "interpretation": proximity_analysis.interpretation
        },

        "properties": {
            "charge_change": property_analysis.charge_change,
            "hydrophobicity_change": property_analysis.hydrophobicity_change,
            "volume_change": property_analysis.volume_change,
            "aromatic_change": property_analysis.aromatic_change,
            "severity_score": property_analysis.severity_score,
            "interpretation": property_analysis.interpretation
        },

        "combined_risk": {
            "score": combined_risk,
            "level": risk_level,
            "interpretation": risk_interpretation
        }
    }


# ============================================================================
# Batch Analysis
# ============================================================================

def analyze_mutations_batch(
    pdb_path: str,
    mutations: List[str],
    chain_start: int = 1,
    sequence_length: int = None,
    is_canonical: bool = True
) -> List[Dict[str, Any]]:
    """
    Analyze multiple mutations against the same structure.

    Args:
        pdb_path: Path to predicted structure PDB
        mutations: List of mutations like ["G324A", "K262R"]
        chain_start: Starting residue of input sequence
        sequence_length: Length of input sequence
        is_canonical: Whether mutations use canonical numbering

    Returns:
        List of analysis results for each mutation
    """
    results = []
    for mutation in mutations:
        try:
            result = analyze_mutation_impact(
                pdb_path, mutation, chain_start, sequence_length, is_canonical
            )
            results.append(result)
        except Exception as e:
            results.append({
                "mutation": mutation,
                "error": str(e)
            })
    return results


# ============================================================================
# Convenience Functions
# ============================================================================

def quick_proximity_check(
    pdb_path: str,
    mutation_pdb_position: int
) -> Dict[str, Any]:
    """
    Quick check of mutation proximity to binding site.

    Args:
        pdb_path: Path to PDB file
        mutation_pdb_position: Residue number in PDB

    Returns:
        Simple dictionary with proximity info
    """
    analysis = analyze_binding_proximity(pdb_path, mutation_pdb_position)

    return {
        "distance": analysis.min_distance_to_ligand,
        "category": analysis.category.value,
        "in_binding_pocket": analysis.is_in_binding_pocket,
        "interpretation": analysis.interpretation
    }


def quick_property_check(wt: str, mut: str) -> Dict[str, Any]:
    """
    Quick check of mutation property changes.

    Args:
        wt: Wild-type residue
        mut: Mutant residue

    Returns:
        Simple dictionary with property changes
    """
    analysis = analyze_mutation_properties(wt, mut)

    return {
        "charge_change": analysis.charge_change,
        "hydrophobicity_change": analysis.hydrophobicity_change,
        "volume_change": analysis.volume_change,
        "severity": analysis.severity_score,
        "interpretation": analysis.interpretation
    }


# ============================================================================
# Main / Demo
# ============================================================================

if __name__ == "__main__":
    import sys

    print("=" * 60)
    print("BoltzOmics Mutation Analysis Module")
    print("=" * 60)

    # Demo property analysis
    print("\nProperty Analysis Examples:")
    for mutation in ["K262R", "K262E", "F656A", "G324A"]:
        wt, pos, mut = parse_mutation(mutation)
        result = analyze_mutation_properties(wt, mut, pos)
        print(f"  {mutation}: {result.interpretation}")

    # Demo with PDB if provided
    if len(sys.argv) > 2:
        pdb_path = sys.argv[1]
        mutation = sys.argv[2]
        chain_start = int(sys.argv[3]) if len(sys.argv) > 3 else 1

        print(f"\nFull Analysis: {mutation} in {pdb_path}")
        result = analyze_mutation_impact(pdb_path, mutation, chain_start)

        print(f"\nProximity: {result['proximity']['interpretation']}")
        print(f"Properties: {result['properties']['interpretation']}")
        print(f"Combined Risk: {result['combined_risk']['level']} ({result['combined_risk']['score']:.2f})")
    else:
        print("\nUsage: python mutation_analysis.py <pdb_file> <mutation> [chain_start]")
        print("Example: python mutation_analysis.py model.pdb K262R 1")

    print("=" * 60)
