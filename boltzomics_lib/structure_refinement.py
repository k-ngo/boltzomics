#!/usr/bin/env python3
"""
Structure Refinement Module for BoltzOmics

This module provides post-processing refinement and analysis of Boltz-2 predicted
structures, including:

1. OpenMM-based energy minimization
2. Binding interface analysis (contacts, clashes, hydrogen bonds)
3. Interaction fingerprint comparison between WT and mutant

These refinements add scientific rigor and provide additional metrics for
interpreting pharmacogenomic predictions.

Author: BoltzOmics Team
"""

import os
import json
import logging
import tempfile
import re
import math
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
from enum import Enum
import numpy as np

# Configure logging
logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes
# ============================================================================

class InteractionType(Enum):
    """Types of protein-ligand interactions."""
    HYDROPHOBIC = "hydrophobic"
    HYDROGEN_BOND = "hydrogen_bond"
    IONIC = "ionic"
    PI_STACKING = "pi_stacking"
    HALOGEN_BOND = "halogen_bond"
    WATER_BRIDGE = "water_bridge"
    CLASH = "clash"


@dataclass
class AtomContact:
    """Represents a contact between protein and ligand atoms."""
    protein_residue: str  # e.g., "ALA123"
    protein_atom: str     # e.g., "CA"
    ligand_atom: str      # e.g., "C1"
    distance: float       # Angstroms
    interaction_type: InteractionType


@dataclass
class BindingInterfaceAnalysis:
    """Results of binding interface analysis."""
    total_contacts: int
    contact_residues: List[str]
    clashes: List[AtomContact]
    hydrogen_bonds: List[AtomContact]
    hydrophobic_contacts: List[AtomContact]
    buried_surface_area: float  # Angstroms^2
    interface_energy_estimate: float  # kcal/mol (approximate)
    clash_score: float  # 0-1, lower is better
    contact_score: float  # 0-1, higher is better


@dataclass
class InteractionFingerprint:
    """Bit vector representing protein-ligand interactions."""
    residue_contacts: Dict[str, bool]  # residue_id -> has_contact
    interaction_types: Dict[str, List[InteractionType]]  # residue_id -> [types]
    fingerprint_bits: List[int]  # Binary fingerprint

    def tanimoto_similarity(self, other: 'InteractionFingerprint') -> float:
        """Calculate Tanimoto similarity between two fingerprints."""
        a = set(i for i, b in enumerate(self.fingerprint_bits) if b)
        b = set(i for i, bit in enumerate(other.fingerprint_bits) if bit)
        if not a and not b:
            return 1.0
        intersection = len(a & b)
        union = len(a | b)
        return intersection / union if union > 0 else 0.0


@dataclass
class RefinementResult:
    """Results of structure refinement."""
    original_pdb_path: str
    refined_pdb_path: Optional[str]
    minimization_performed: bool
    initial_energy: Optional[float]  # kcal/mol
    final_energy: Optional[float]    # kcal/mol
    energy_change: Optional[float]   # kcal/mol
    rmsd_from_original: Optional[float]  # Angstroms
    clashes_removed: int
    refinement_time_seconds: float


@dataclass
class MutationInterfaceComparison:
    """Comparison of binding interfaces between WT and mutant."""
    wt_analysis: BindingInterfaceAnalysis
    mutant_analysis: BindingInterfaceAnalysis
    wt_fingerprint: InteractionFingerprint
    mutant_fingerprint: InteractionFingerprint
    fingerprint_similarity: float  # Tanimoto similarity
    lost_contacts: List[str]  # Residues that lost contact
    gained_contacts: List[str]  # Residues that gained contact
    new_clashes: List[AtomContact]  # Clashes introduced by mutation
    delta_buried_area: float  # Change in buried surface area
    overall_disruption_score: float  # 0-1, higher = more disruption


# ============================================================================
# OpenMM Energy Minimization
# ============================================================================

def check_openmm_available() -> bool:
    """Check if OpenMM is installed and available."""
    try:
        import openmm
        import openmm.app
        return True
    except ImportError:
        return False


def minimize_structure_openmm(
    pdb_path: str,
    output_path: Optional[str] = None,
    max_iterations: int = 500,
    tolerance: float = 10.0,  # kJ/mol/nm
    add_hydrogens: bool = True,
    fix_residues: bool = True,
    platform: str = "CPU"
) -> RefinementResult:
    """
    Perform energy minimization on a PDB structure using OpenMM.

    This relaxes the structure to remove clashes and improve geometry
    while preserving the overall fold predicted by Boltz-2.

    Args:
        pdb_path: Path to input PDB file
        output_path: Path for refined PDB (default: input_minimized.pdb)
        max_iterations: Maximum minimization steps
        tolerance: Energy tolerance for convergence
        add_hydrogens: Whether to add missing hydrogens
        fix_residues: Whether to attempt to fix non-standard residues
        platform: OpenMM platform (CPU, CUDA, OpenCL)

    Returns:
        RefinementResult with minimization details
    """
    import time
    start_time = time.time()

    if not check_openmm_available():
        logger.warning("OpenMM not available. Skipping energy minimization.")
        return RefinementResult(
            original_pdb_path=pdb_path,
            refined_pdb_path=None,
            minimization_performed=False,
            initial_energy=None,
            final_energy=None,
            energy_change=None,
            rmsd_from_original=None,
            clashes_removed=0,
            refinement_time_seconds=time.time() - start_time
        )

    try:
        from openmm import app, unit, LangevinMiddleIntegrator, Platform
        from openmm.app import PDBFile, ForceField, Modeller, PDBxFile
        import openmm

        # Set output path
        if output_path is None:
            base = os.path.splitext(pdb_path)[0]
            output_path = f"{base}_minimized.pdb"

        # Load the structure
        logger.info(f"Loading structure from {pdb_path}")
        pdb = PDBFile(pdb_path)

        # Create modeller for preprocessing
        modeller = Modeller(pdb.topology, pdb.positions)
        removed_residues: List[Tuple[int, str, str, str]] = []

        # Use Amber force field (good for proteins)
        # For ligands, we need GAFF - use implicit solvent for simplicity
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

        # Add hydrogens if requested
        if add_hydrogens:
            logger.info("Adding missing hydrogens")
            # Use OpenMM's built-in hydrogen placement first to avoid immediate
            # template matching failures on partially nonstandard topologies.
            modeller.addHydrogens()

        # Create system with implicit solvent (faster than explicit)
        logger.info("Creating OpenMM system")
        system = None
        max_template_retries = 24
        for _ in range(max_template_retries):
            try:
                system = forcefield.createSystem(
                    modeller.topology,
                    nonbondedMethod=app.NoCutoff,
                    constraints=app.HBonds,
                    ignoreExternalBonds=True
                )
                break
            except Exception as create_exc:
                template_idx = _extract_template_error_residue_index(str(create_exc))
                residues = list(modeller.topology.residues())
                if template_idx is None or template_idx < 0 or template_idx >= len(residues):
                    raise
                bad_res = residues[template_idx]
                removed_residues.append((template_idx, bad_res.name, bad_res.id, bad_res.chain.id))
                logger.warning(
                    "Removing unparameterized residue during minimization fallback: "
                    "index=%s name=%s id=%s chain=%s",
                    template_idx,
                    bad_res.name,
                    bad_res.id,
                    bad_res.chain.id,
                )
                modeller.delete([bad_res])

        if system is None:
            raise RuntimeError("Failed to create OpenMM system after template fallback retries.")

        # Create integrator (not used for minimization but required)
        integrator = LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,
            0.002 * unit.picoseconds
        )

        # Select platform
        try:
            platform_obj = Platform.getPlatformByName(platform)
        except Exception:
            platform_obj = Platform.getPlatformByName("CPU")

        # Create simulation
        simulation = app.Simulation(
            modeller.topology,
            system,
            integrator,
            platform_obj
        )
        simulation.context.setPositions(modeller.positions)

        # Get initial energy
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        initial_positions = state.getPositions(asNumpy=True)
        final_energy = initial_energy
        final_positions = initial_positions

        # Optional ligand-aware safety guard:
        # run minimization in chunks and stop/revert if protein-ligand clashes worsen.
        ligand_positions_angstrom = _extract_ligand_positions_from_pdb(pdb_path)
        use_ligand_guard = ligand_positions_angstrom is not None and ligand_positions_angstrom.shape[0] > 0
        baseline_ligand_metrics = None
        if use_ligand_guard:
            baseline_ligand_metrics = _protein_ligand_clash_metrics(
                positions=initial_positions,
                topology=modeller.topology,
                ligand_positions_angstrom=ligand_positions_angstrom,
                unit_module=unit,
            )

        # Count initial clashes
        initial_clashes = _count_clashes(initial_positions, modeller.topology)

        # Perform minimization (chunked mode with ligand clash guard when ligand exists).
        logger.info(f"Running energy minimization (max {max_iterations} iterations)")
        if use_ligand_guard and int(max_iterations) > 0:
            chunk_count = 5
            chunk_iters = max(1, int(math.ceil(float(max_iterations) / float(chunk_count))))
            remaining = int(max_iterations)
            last_safe_positions = initial_positions
            last_safe_energy = initial_energy
            guard_stopped = False

            while remaining > 0:
                step_iters = min(chunk_iters, remaining)
                simulation.minimizeEnergy(
                    tolerance=tolerance * unit.kilojoules_per_mole / unit.nanometer,
                    maxIterations=step_iters,
                )
                state = simulation.context.getState(getEnergy=True, getPositions=True)
                candidate_positions = state.getPositions(asNumpy=True)
                candidate_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

                candidate_metrics = _protein_ligand_clash_metrics(
                    positions=candidate_positions,
                    topology=modeller.topology,
                    ligand_positions_angstrom=ligand_positions_angstrom,
                    unit_module=unit,
                )

                # Define worsening by strictly increased severe-close-contact counts.
                baseline_n_lt_2A = int((baseline_ligand_metrics or {}).get("n_lt_2A", 0))
                baseline_n_lt_1_6A = int((baseline_ligand_metrics or {}).get("n_lt_1_6A", 0))
                candidate_n_lt_2A = int(candidate_metrics.get("n_lt_2A", 0))
                candidate_n_lt_1_6A = int(candidate_metrics.get("n_lt_1_6A", 0))
                worsened = (
                    candidate_n_lt_1_6A > baseline_n_lt_1_6A
                    or candidate_n_lt_2A > baseline_n_lt_2A
                )

                if worsened:
                    logger.warning(
                        "Stopping minimization early due to worsened protein-ligand clashes: "
                        "baseline(<1.6A=%d,<2.0A=%d) -> candidate(<1.6A=%d,<2.0A=%d)",
                        baseline_n_lt_1_6A,
                        baseline_n_lt_2A,
                        candidate_n_lt_1_6A,
                        candidate_n_lt_2A,
                    )
                    simulation.context.setPositions(last_safe_positions)
                    guard_stopped = True
                    break

                last_safe_positions = candidate_positions
                last_safe_energy = candidate_energy
                remaining -= step_iters

            final_positions = last_safe_positions
            final_energy = last_safe_energy
            if guard_stopped:
                logger.info("Ligand clash guard applied: reverted to last safe minimization state.")
        else:
            simulation.minimizeEnergy(
                tolerance=tolerance * unit.kilojoules_per_mole / unit.nanometer,
                maxIterations=max_iterations
            )
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            final_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            final_positions = state.getPositions(asNumpy=True)

        # Count final clashes
        final_clashes = _count_clashes(final_positions, modeller.topology)

        # Calculate RMSD
        rmsd = _calculate_rmsd(initial_positions, final_positions)

        # Save minimized structure
        logger.info(f"Saving minimized structure to {output_path}")
        if removed_residues:
            _write_minimized_coords_preserve_original_records(
                original_pdb_path=pdb_path,
                output_pdb_path=output_path,
                topology=modeller.topology,
                positions=final_positions,
                unit_module=unit,
            )
            logger.info(
                "Preserved original records while writing minimized coordinates "
                "after removing %d unparameterized residues for minimization.",
                len(removed_residues),
            )
        else:
            with open(output_path, 'w') as f:
                PDBFile.writeFile(modeller.topology, final_positions, f)

        return RefinementResult(
            original_pdb_path=pdb_path,
            refined_pdb_path=output_path,
            minimization_performed=True,
            initial_energy=initial_energy,
            final_energy=final_energy,
            energy_change=final_energy - initial_energy,
            rmsd_from_original=rmsd,
            clashes_removed=max(0, initial_clashes - final_clashes),
            refinement_time_seconds=time.time() - start_time
        )

    except Exception as e:
        logger.error(f"OpenMM minimization failed: {e}")
        return RefinementResult(
            original_pdb_path=pdb_path,
            refined_pdb_path=None,
            minimization_performed=False,
            initial_energy=None,
            final_energy=None,
            energy_change=None,
            rmsd_from_original=None,
            clashes_removed=0,
            refinement_time_seconds=time.time() - start_time
        )


def _extract_template_error_residue_index(error_text: str) -> Optional[int]:
    """
    Parse OpenMM template matching errors and return residue index if available.

    Expected pattern:
    "No template found for residue <idx> (RES)..."
    """
    match = re.search(r"No template found for residue\s+(\d+)\s+\(", error_text)
    if not match:
        return None
    try:
        return int(match.group(1))
    except Exception:
        return None


def _make_atom_coord_key(atom) -> Tuple[str, str, str, str]:
    """Build a stable coordinate key from OpenMM topology atom metadata."""
    residue = atom.residue
    chain_id = str(residue.chain.id)
    residue_id = str(residue.id).strip()
    residue_name = str(residue.name).strip()
    atom_name = str(atom.name).strip()
    return (chain_id, residue_id, residue_name, atom_name)


def _write_minimized_coords_preserve_original_records(
    original_pdb_path: str,
    output_pdb_path: str,
    topology,
    positions,
    unit_module,
) -> None:
    """
    Write a minimized PDB while preserving original records not parameterized by OpenMM.

    This keeps unsupported residues (for example ligands) in the output unchanged,
    while updating coordinates for atoms included in minimization.
    """
    # Positions are in nm in OpenMM; convert to Angstrom for PDB formatting.
    positions_angstrom = positions.value_in_unit(unit_module.angstrom)
    coord_map: Dict[Tuple[str, str, str, str], Tuple[float, float, float]] = {}
    for atom, xyz in zip(topology.atoms(), positions_angstrom):
        key = _make_atom_coord_key(atom)
        coord_map[key] = (float(xyz[0]), float(xyz[1]), float(xyz[2]))

    with open(original_pdb_path, "r", encoding="utf-8", errors="ignore") as src, \
            open(output_pdb_path, "w", encoding="utf-8") as dst:
        for line in src:
            record = line[:6]
            if record in ("ATOM  ", "HETATM"):
                if len(line) >= 54:
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain_id = (line[21] or " ").strip()
                    residue_id = line[22:26].strip()
                    key = (chain_id, residue_id, residue_name, atom_name)
                    coords = coord_map.get(key)
                    if coords is not None:
                        x, y, z = coords
                        line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            dst.write(line)


def _count_clashes(positions: np.ndarray, topology, clash_threshold: float = 0.4) -> int:
    """Count atomic clashes (distances below threshold in nm)."""
    n_atoms = len(positions)
    clash_count = 0

    # Get positions as numpy array
    if hasattr(positions, 'value_in_unit'):
        from openmm import unit
        pos = positions.value_in_unit(unit.nanometer)
    else:
        pos = positions

    # Simple O(n^2) clash detection - could be optimized with spatial hashing
    for i in range(min(n_atoms, 1000)):  # Limit to first 1000 atoms for speed
        for j in range(i + 4, min(n_atoms, 1000)):  # Skip bonded atoms
            dist = np.linalg.norm(pos[i] - pos[j])
            if dist < clash_threshold:
                clash_count += 1

    return clash_count


def _extract_ligand_positions_from_pdb(pdb_path: str, ligand_resname: str = "LIG") -> Optional[np.ndarray]:
    """Extract ligand atom coordinates (Angstrom) from a PDB file."""
    coords: List[Tuple[float, float, float]] = []
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("HETATM"):
                    continue
                if line[17:20].strip() != ligand_resname:
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    continue
                coords.append((x, y, z))
    except Exception:
        return None
    if not coords:
        return None
    return np.asarray(coords, dtype=float)


def _protein_ligand_clash_metrics(
    positions,
    topology,
    ligand_positions_angstrom: np.ndarray,
    unit_module,
) -> Dict[str, float]:
    """
    Compute protein-ligand proximity metrics from current simulation positions.

    Returns:
        min distance and close-contact counts at <2.0A and <1.6A.
    """
    if hasattr(positions, "value_in_unit"):
        positions_angstrom = positions.value_in_unit(unit_module.angstrom)
    else:
        positions_angstrom = positions

    protein_coords: List[Tuple[float, float, float]] = []
    for atom in topology.atoms():
        if atom.element is None:
            continue
        if atom.element.symbol == "H":
            continue
        if atom.residue.name == "LIG":
            continue
        xyz = positions_angstrom[atom.index]
        protein_coords.append((float(xyz[0]), float(xyz[1]), float(xyz[2])))

    if not protein_coords or ligand_positions_angstrom is None or ligand_positions_angstrom.size == 0:
        return {"min_dist_A": float("inf"), "n_lt_2A": 0, "n_lt_1_6A": 0}

    prot = np.asarray(protein_coords, dtype=float)
    lig = np.asarray(ligand_positions_angstrom, dtype=float)
    dists = np.sqrt(((prot[:, None, :] - lig[None, :, :]) ** 2).sum(axis=2))
    return {
        "min_dist_A": float(np.min(dists)),
        "n_lt_2A": int(np.sum(dists < 2.0)),
        "n_lt_1_6A": int(np.sum(dists < 1.6)),
    }


def _calculate_rmsd(pos1: np.ndarray, pos2: np.ndarray) -> float:
    """Calculate RMSD between two position arrays."""
    if hasattr(pos1, 'value_in_unit'):
        from openmm import unit
        pos1 = pos1.value_in_unit(unit.angstrom)
        pos2 = pos2.value_in_unit(unit.angstrom)

    diff = pos1 - pos2
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


# ============================================================================
# Binding Interface Analysis
# ============================================================================

def analyze_binding_interface(
    pdb_path: str,
    ligand_chain: str = "B",
    protein_chain: str = "A",
    contact_threshold: float = 4.5,  # Angstroms
    clash_threshold: float = 2.0,    # Angstroms
    hbond_threshold: float = 3.5     # Angstroms
) -> BindingInterfaceAnalysis:
    """
    Analyze the protein-ligand binding interface.

    Identifies contacts, clashes, hydrogen bonds, and estimates
    interface properties.

    Args:
        pdb_path: Path to PDB file
        ligand_chain: Chain ID for ligand
        protein_chain: Chain ID for protein
        contact_threshold: Distance cutoff for contacts
        clash_threshold: Distance cutoff for clashes
        hbond_threshold: Distance cutoff for H-bonds

    Returns:
        BindingInterfaceAnalysis with detailed interface metrics
    """
    # Parse PDB file
    protein_atoms, ligand_atoms = _parse_pdb_atoms(pdb_path, protein_chain, ligand_chain)

    if not protein_atoms or not ligand_atoms:
        logger.warning(f"Could not find protein ({protein_chain}) or ligand ({ligand_chain}) in {pdb_path}")
        return BindingInterfaceAnalysis(
            total_contacts=0,
            contact_residues=[],
            clashes=[],
            hydrogen_bonds=[],
            hydrophobic_contacts=[],
            buried_surface_area=0.0,
            interface_energy_estimate=0.0,
            clash_score=0.0,
            contact_score=0.0
        )

    # Find all contacts
    contacts = []
    clashes = []
    hydrogen_bonds = []
    hydrophobic_contacts = []
    contact_residues = set()

    # H-bond donors/acceptors
    hbond_atoms = {'N', 'O', 'NE', 'NE1', 'NE2', 'ND1', 'ND2', 'NZ', 'OG', 'OG1', 'OH', 'OD1', 'OD2', 'OE1', 'OE2'}
    hydrophobic_atoms = {'C', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CZ'}

    for p_atom in protein_atoms:
        for l_atom in ligand_atoms:
            dist = np.linalg.norm(np.array(p_atom['coord']) - np.array(l_atom['coord']))

            if dist < contact_threshold:
                residue_id = f"{p_atom['resname']}{p_atom['resid']}"
                contact_residues.add(residue_id)

                # Determine interaction type
                if dist < clash_threshold:
                    interaction_type = InteractionType.CLASH
                    contact = AtomContact(
                        protein_residue=residue_id,
                        protein_atom=p_atom['name'],
                        ligand_atom=l_atom['name'],
                        distance=dist,
                        interaction_type=interaction_type
                    )
                    clashes.append(contact)
                    contacts.append(contact)

                elif p_atom['name'] in hbond_atoms and dist < hbond_threshold:
                    interaction_type = InteractionType.HYDROGEN_BOND
                    contact = AtomContact(
                        protein_residue=residue_id,
                        protein_atom=p_atom['name'],
                        ligand_atom=l_atom['name'],
                        distance=dist,
                        interaction_type=interaction_type
                    )
                    hydrogen_bonds.append(contact)
                    contacts.append(contact)

                elif p_atom['name'] in hydrophobic_atoms:
                    interaction_type = InteractionType.HYDROPHOBIC
                    contact = AtomContact(
                        protein_residue=residue_id,
                        protein_atom=p_atom['name'],
                        ligand_atom=l_atom['name'],
                        distance=dist,
                        interaction_type=interaction_type
                    )
                    hydrophobic_contacts.append(contact)
                    contacts.append(contact)

                else:
                    contact = AtomContact(
                        protein_residue=residue_id,
                        protein_atom=p_atom['name'],
                        ligand_atom=l_atom['name'],
                        distance=dist,
                        interaction_type=InteractionType.HYDROPHOBIC  # Default
                    )
                    contacts.append(contact)

    # Estimate buried surface area (simplified)
    bsa = len(contact_residues) * 120.0  # ~120 A^2 per contacting residue (rough estimate)

    # Estimate interface energy (very rough)
    energy_estimate = -0.5 * len(hydrogen_bonds) - 0.1 * len(hydrophobic_contacts) + 2.0 * len(clashes)

    # Calculate scores
    clash_score = len(clashes) / max(len(contacts), 1)
    contact_score = (len(hydrogen_bonds) * 2 + len(hydrophobic_contacts)) / max(len(contacts), 1)

    return BindingInterfaceAnalysis(
        total_contacts=len(contacts),
        contact_residues=list(contact_residues),
        clashes=clashes,
        hydrogen_bonds=hydrogen_bonds,
        hydrophobic_contacts=hydrophobic_contacts,
        buried_surface_area=bsa,
        interface_energy_estimate=energy_estimate,
        clash_score=min(1.0, clash_score),
        contact_score=min(1.0, contact_score)
    )


def detect_chains(pdb_path: str) -> Tuple[str, str]:
    """
    Auto-detect protein and ligand chain IDs from PDB file.

    Assumes protein is the largest chain (by atom count) using ATOM records,
    and ligand is HETATM or a smaller chain.

    Returns:
        Tuple of (protein_chain, ligand_chain)
    """
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
                # Skip water and common ions
                if resname not in ['HOH', 'WAT', 'NA', 'CL', 'MG', 'ZN', 'CA']:
                    chain_hetatm_counts[chain] = chain_hetatm_counts.get(chain, 0) + 1

    # Protein is the largest ATOM chain
    protein_chain = max(chain_atom_counts, key=chain_atom_counts.get) if chain_atom_counts else 'A'

    # Ligand is preferentially HETATM, or the next largest chain
    if chain_hetatm_counts:
        # Prefer HETATM chain that's not the protein
        non_protein_hetatm = {k: v for k, v in chain_hetatm_counts.items() if k != protein_chain}
        if non_protein_hetatm:
            ligand_chain = max(non_protein_hetatm, key=non_protein_hetatm.get)
        else:
            ligand_chain = max(chain_hetatm_counts, key=chain_hetatm_counts.get)
    else:
        # Fall back to second largest ATOM chain
        other_chains = {k: v for k, v in chain_atom_counts.items() if k != protein_chain}
        ligand_chain = max(other_chains, key=other_chains.get) if other_chains else 'B'

    return protein_chain, ligand_chain


def _parse_pdb_atoms(
    pdb_path: str,
    protein_chain: str,
    ligand_chain: str
) -> Tuple[List[Dict], List[Dict]]:
    """Parse atoms from PDB file."""
    protein_atoms = []
    ligand_atoms = []

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
                    protein_atoms.append(atom)
                elif atom['chain'] == ligand_chain:
                    ligand_atoms.append(atom)

    return protein_atoms, ligand_atoms


# ============================================================================
# Interaction Fingerprint
# ============================================================================

def generate_interaction_fingerprint(
    pdb_path: str,
    ligand_chain: str = "B",
    protein_chain: str = "A",
    contact_threshold: float = 4.5,
    max_residues: int = 500
) -> InteractionFingerprint:
    """
    Generate an interaction fingerprint for a protein-ligand complex.

    The fingerprint encodes which residues are in contact with the ligand
    and what types of interactions they make.

    Args:
        pdb_path: Path to PDB file
        ligand_chain: Chain ID for ligand
        protein_chain: Chain ID for protein
        contact_threshold: Distance cutoff for contacts
        max_residues: Maximum residue number to consider

    Returns:
        InteractionFingerprint with contact information
    """
    # Get interface analysis
    analysis = analyze_binding_interface(
        pdb_path, ligand_chain, protein_chain, contact_threshold
    )

    # Build residue contact map
    residue_contacts = {}
    interaction_types = {}

    for contact in (analysis.hydrogen_bonds + analysis.hydrophobic_contacts + analysis.clashes):
        res = contact.protein_residue
        residue_contacts[res] = True
        if res not in interaction_types:
            interaction_types[res] = []
        if contact.interaction_type not in interaction_types[res]:
            interaction_types[res].append(contact.interaction_type)

    # Generate binary fingerprint
    # Format: [residue_has_contact] + [residue_has_hbond] + [residue_has_hydrophobic] + [residue_has_clash]
    fingerprint_bits = []

    for i in range(1, max_residues + 1):
        # Find residue by number
        res_match = None
        for res in residue_contacts:
            try:
                if int(''.join(filter(str.isdigit, res))) == i:
                    res_match = res
                    break
            except ValueError:
                continue

        # Contact bit
        has_contact = res_match is not None
        fingerprint_bits.append(1 if has_contact else 0)

        if has_contact and res_match in interaction_types:
            types = interaction_types[res_match]
            fingerprint_bits.append(1 if InteractionType.HYDROGEN_BOND in types else 0)
            fingerprint_bits.append(1 if InteractionType.HYDROPHOBIC in types else 0)
            fingerprint_bits.append(1 if InteractionType.CLASH in types else 0)
        else:
            fingerprint_bits.extend([0, 0, 0])

    return InteractionFingerprint(
        residue_contacts=residue_contacts,
        interaction_types=interaction_types,
        fingerprint_bits=fingerprint_bits
    )


# ============================================================================
# WT vs Mutant Comparison
# ============================================================================

def compare_mutation_interfaces(
    wt_pdb_path: str,
    mutant_pdb_path: str,
    mutation_position: int,
    ligand_chain: str = "B",
    protein_chain: str = "A"
) -> MutationInterfaceComparison:
    """
    Compare binding interfaces between wild-type and mutant structures.

    This analysis identifies how a mutation affects:
    - Contact residues (lost/gained contacts)
    - Interaction types (H-bonds, hydrophobic)
    - Steric clashes (introduced by mutation)
    - Overall interface similarity

    Args:
        wt_pdb_path: Path to wild-type PDB
        mutant_pdb_path: Path to mutant PDB
        mutation_position: Residue number of the mutation
        ligand_chain: Chain ID for ligand
        protein_chain: Chain ID for protein

    Returns:
        MutationInterfaceComparison with detailed comparison
    """
    # Analyze both structures
    wt_analysis = analyze_binding_interface(wt_pdb_path, ligand_chain, protein_chain)
    mut_analysis = analyze_binding_interface(mutant_pdb_path, ligand_chain, protein_chain)

    # Generate fingerprints
    wt_fp = generate_interaction_fingerprint(wt_pdb_path, ligand_chain, protein_chain)
    mut_fp = generate_interaction_fingerprint(mutant_pdb_path, ligand_chain, protein_chain)

    # Calculate fingerprint similarity
    fp_similarity = wt_fp.tanimoto_similarity(mut_fp)

    # Find lost and gained contacts
    wt_residues = set(wt_analysis.contact_residues)
    mut_residues = set(mut_analysis.contact_residues)

    lost_contacts = list(wt_residues - mut_residues)
    gained_contacts = list(mut_residues - wt_residues)

    # Identify new clashes
    wt_clash_residues = {c.protein_residue for c in wt_analysis.clashes}
    new_clashes = [c for c in mut_analysis.clashes if c.protein_residue not in wt_clash_residues]

    # Calculate overall disruption score
    # Higher = more disruption
    disruption = 0.0

    # Contact changes
    total_contacts = len(wt_residues | mut_residues)
    if total_contacts > 0:
        contact_change = len(lost_contacts) + len(gained_contacts)
        disruption += 0.3 * (contact_change / total_contacts)

    # New clashes are highly disruptive
    disruption += 0.4 * min(1.0, len(new_clashes) / 3.0)

    # Fingerprint dissimilarity
    disruption += 0.3 * (1.0 - fp_similarity)

    # Check if mutation position is in contact
    mut_pos_str = str(mutation_position)
    mut_in_contact_wt = any(mut_pos_str in res for res in wt_residues)
    mut_in_contact_mut = any(mut_pos_str in res for res in mut_residues)

    if mut_in_contact_wt or mut_in_contact_mut:
        # Mutation at interface is more impactful
        disruption = min(1.0, disruption * 1.5)

    return MutationInterfaceComparison(
        wt_analysis=wt_analysis,
        mutant_analysis=mut_analysis,
        wt_fingerprint=wt_fp,
        mutant_fingerprint=mut_fp,
        fingerprint_similarity=fp_similarity,
        lost_contacts=lost_contacts,
        gained_contacts=gained_contacts,
        new_clashes=new_clashes,
        delta_buried_area=mut_analysis.buried_surface_area - wt_analysis.buried_surface_area,
        overall_disruption_score=min(1.0, disruption)
    )


# ============================================================================
# RDKit-based Enhanced Analysis
# ============================================================================

def analyze_ligand_interactions_rdkit(
    pdb_path: str,
    ligand_smiles: str,
    ligand_chain: str = "B"
) -> Dict[str, Any]:
    """
    Enhanced interaction analysis using RDKit.

    Provides more sophisticated interaction detection including:
    - Pi-stacking
    - Halogen bonds
    - Salt bridges

    Args:
        pdb_path: Path to PDB file
        ligand_smiles: SMILES string of the ligand
        ligand_chain: Chain ID for ligand

    Returns:
        Dictionary with detailed interaction information
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError:
        logger.warning("RDKit not available for enhanced analysis")
        return {"error": "RDKit not available"}

    # Load ligand from SMILES
    mol = Chem.MolFromSmiles(ligand_smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}

    # Add hydrogens and generate 3D conformer
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    # Calculate molecular properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rotatable = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)

    # Count aromatic rings (for pi-stacking potential)
    aromatic_rings = sum(1 for ring in mol.GetRingInfo().AtomRings()
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring))

    # Check for halogens (for halogen bonding)
    halogens = sum(1 for atom in mol.GetAtoms()
                   if atom.GetAtomicNum() in [9, 17, 35, 53])  # F, Cl, Br, I

    return {
        "molecular_weight": mw,
        "logp": logp,
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": rotatable,
        "tpsa": tpsa,
        "aromatic_rings": aromatic_rings,
        "halogens": halogens,
        "pi_stacking_potential": aromatic_rings > 0,
        "halogen_bonding_potential": halogens > 0
    }


# ============================================================================
# Combined Refinement Pipeline
# ============================================================================

def refine_and_analyze(
    pdb_path: str,
    ligand_smiles: Optional[str] = None,
    ligand_chain: str = "B",
    protein_chain: str = "A",
    perform_minimization: bool = True,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    Complete refinement and analysis pipeline.

    Performs:
    1. Optional OpenMM energy minimization
    2. Binding interface analysis
    3. Interaction fingerprint generation
    4. Optional RDKit-based ligand analysis

    Args:
        pdb_path: Path to input PDB
        ligand_smiles: Optional SMILES for enhanced analysis
        ligand_chain: Chain ID for ligand
        protein_chain: Chain ID for protein
        perform_minimization: Whether to minimize structure
        output_dir: Directory for output files

    Returns:
        Dictionary with all analysis results
    """
    results = {
        "input_pdb": pdb_path,
        "success": True
    }

    # Set up output directory
    if output_dir is None:
        output_dir = os.path.dirname(pdb_path)

    # 1. Optional minimization
    if perform_minimization and check_openmm_available():
        logger.info("Performing energy minimization...")
        output_pdb = os.path.join(
            output_dir,
            os.path.basename(pdb_path).replace('.pdb', '_minimized.pdb')
        )
        refinement = minimize_structure_openmm(pdb_path, output_pdb)
        results["refinement"] = {
            "performed": refinement.minimization_performed,
            "refined_pdb": refinement.refined_pdb_path,
            "initial_energy": refinement.initial_energy,
            "final_energy": refinement.final_energy,
            "energy_change": refinement.energy_change,
            "rmsd": refinement.rmsd_from_original,
            "clashes_removed": refinement.clashes_removed,
            "time_seconds": refinement.refinement_time_seconds
        }
        # Use refined structure for analysis if available
        analysis_pdb = refinement.refined_pdb_path or pdb_path
    else:
        results["refinement"] = {"performed": False, "reason": "minimization disabled or OpenMM unavailable"}
        analysis_pdb = pdb_path

    # 2. Interface analysis
    logger.info("Analyzing binding interface...")
    interface = analyze_binding_interface(analysis_pdb, ligand_chain, protein_chain)
    results["interface"] = {
        "total_contacts": interface.total_contacts,
        "contact_residues": interface.contact_residues,
        "n_clashes": len(interface.clashes),
        "n_hydrogen_bonds": len(interface.hydrogen_bonds),
        "n_hydrophobic_contacts": len(interface.hydrophobic_contacts),
        "buried_surface_area": interface.buried_surface_area,
        "interface_energy_estimate": interface.interface_energy_estimate,
        "clash_score": interface.clash_score,
        "contact_score": interface.contact_score
    }

    # 3. Interaction fingerprint
    logger.info("Generating interaction fingerprint...")
    fingerprint = generate_interaction_fingerprint(analysis_pdb, ligand_chain, protein_chain)
    results["fingerprint"] = {
        "contact_residues": list(fingerprint.residue_contacts.keys()),
        "n_contact_residues": len(fingerprint.residue_contacts),
        "fingerprint_length": len(fingerprint.fingerprint_bits),
        "fingerprint_density": sum(fingerprint.fingerprint_bits) / len(fingerprint.fingerprint_bits) if fingerprint.fingerprint_bits else 0
    }

    # 4. RDKit analysis if SMILES provided
    if ligand_smiles:
        logger.info("Performing RDKit ligand analysis...")
        rdkit_analysis = analyze_ligand_interactions_rdkit(pdb_path, ligand_smiles, ligand_chain)
        results["ligand_properties"] = rdkit_analysis

    # Save results to JSON
    results_path = os.path.join(
        output_dir,
        os.path.basename(pdb_path).replace('.pdb', '_refinement_analysis.json')
    )
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    results["results_file"] = results_path

    logger.info(f"Refinement and analysis complete. Results saved to {results_path}")
    return results


# ============================================================================
# Batch Processing
# ============================================================================

def batch_refine_screening_results(
    results_dir: str,
    perform_minimization: bool = False,  # Default off for speed
    ligand_chain: str = "B",
    protein_chain: str = "A"
) -> List[Dict[str, Any]]:
    """
    Batch process all PDB files from a screening directory.

    Args:
        results_dir: Directory containing screening results
        perform_minimization: Whether to perform energy minimization
        ligand_chain: Chain ID for ligand
        protein_chain: Chain ID for protein

    Returns:
        List of analysis results for each structure
    """
    from pathlib import Path

    results = []
    pdb_files = list(Path(results_dir).rglob("*.pdb"))

    logger.info(f"Found {len(pdb_files)} PDB files to process")

    for pdb_path in pdb_files:
        if "_minimized" in str(pdb_path):
            continue  # Skip already minimized files

        logger.info(f"Processing {pdb_path}")
        try:
            result = refine_and_analyze(
                str(pdb_path),
                ligand_chain=ligand_chain,
                protein_chain=protein_chain,
                perform_minimization=perform_minimization
            )
            results.append(result)
        except Exception as e:
            logger.error(f"Failed to process {pdb_path}: {e}")
            results.append({
                "input_pdb": str(pdb_path),
                "success": False,
                "error": str(e)
            })

    return results


# ============================================================================
# Convenience Functions
# ============================================================================

def quick_interface_check(pdb_path: str, auto_detect_chains: bool = True) -> Dict[str, Any]:
    """
    Quick check of binding interface quality.

    Returns a simple summary of interface metrics.
    """
    # Auto-detect chains if requested
    if auto_detect_chains:
        protein_chain, ligand_chain = detect_chains(pdb_path)
    else:
        protein_chain, ligand_chain = "A", "B"

    analysis = analyze_binding_interface(pdb_path, ligand_chain, protein_chain)

    quality = "GOOD"
    if analysis.clash_score > 0.2:
        quality = "POOR - has clashes"
    elif analysis.total_contacts < 5:
        quality = "WEAK - few contacts"
    elif len(analysis.hydrogen_bonds) == 0:
        quality = "MODERATE - no H-bonds"

    return {
        "quality": quality,
        "contacts": analysis.total_contacts,
        "clashes": len(analysis.clashes),
        "h_bonds": len(analysis.hydrogen_bonds),
        "hydrophobic": len(analysis.hydrophobic_contacts),
        "contact_residues": analysis.contact_residues
    }


def compare_wt_mutant_quick(
    wt_pdb: str,
    mutant_pdb: str,
    mutation_position: int,
    auto_detect_chains: bool = True
) -> Dict[str, Any]:
    """
    Quick comparison of WT vs mutant binding.

    Returns summary of how mutation affects binding.
    """
    # Auto-detect chains from WT structure
    if auto_detect_chains:
        protein_chain, ligand_chain = detect_chains(wt_pdb)
    else:
        protein_chain, ligand_chain = "A", "B"

    comparison = compare_mutation_interfaces(
        wt_pdb, mutant_pdb, mutation_position, ligand_chain, protein_chain
    )

    impact = "MINIMAL"
    if comparison.overall_disruption_score > 0.5:
        impact = "MAJOR"
    elif comparison.overall_disruption_score > 0.2:
        impact = "MODERATE"
    elif len(comparison.new_clashes) > 0:
        impact = "MODERATE - introduces clashes"

    return {
        "impact": impact,
        "disruption_score": comparison.overall_disruption_score,
        "fingerprint_similarity": comparison.fingerprint_similarity,
        "lost_contacts": comparison.lost_contacts,
        "gained_contacts": comparison.gained_contacts,
        "new_clashes": len(comparison.new_clashes),
        "delta_buried_area": comparison.delta_buried_area
    }


# ============================================================================
# Main / Demo
# ============================================================================

if __name__ == "__main__":
    import sys

    print("=" * 60)
    print("BoltzOmics Structure Refinement Module")
    print("=" * 60)

    # Check dependencies
    print("\nDependency check:")
    print(f"  OpenMM available: {check_openmm_available()}")

    try:
        from rdkit import Chem
        print("  RDKit available: True")
    except ImportError:
        print("  RDKit available: False")

    # Demo with example if provided
    if len(sys.argv) > 1:
        pdb_path = sys.argv[1]
        print(f"\nAnalyzing: {pdb_path}")

        result = quick_interface_check(pdb_path)
        print(f"\nInterface Quality: {result['quality']}")
        print(f"  Total contacts: {result['contacts']}")
        print(f"  Clashes: {result['clashes']}")
        print(f"  H-bonds: {result['h_bonds']}")
        print(f"  Hydrophobic: {result['hydrophobic']}")
        print(f"  Contact residues: {', '.join(result['contact_residues'][:10])}...")
    else:
        print("\nUsage: python structure_refinement.py <pdb_file>")
        print("\nExample:")
        print("  python structure_refinement.py model_0.pdb")

    print("=" * 60)
