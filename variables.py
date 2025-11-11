"""
Variable descriptions for drug screening visualization metrics.

This module provides human-readable descriptions for the various metrics
used in the CATDiscovery drug screening analysis.
"""

METRIC_DESCRIPTIONS = {
    "pic50": (
        "pIC50 is the negative logarithm (base 10) of the IC50 value. "
        "It provides a more intuitive scale where higher values indicate stronger binding affinity. "
        "For example, an IC50 of 1 μM corresponds to a pIC50 of 6, while an IC50 of 0.001 μM (1 nM) "
        "corresponds to a pIC50 of 9. Typically, pIC50 > 6 indicates good binding affinity."
    ),
    "predicted_ic50": (
        "IC50 (Half-maximal Inhibitory Concentration) is the concentration of a drug required to "
        "inhibit 50% of the target protein's activity. It's measured in micromolar (μM) units. "
        "Lower IC50 values indicate stronger binding and better drug efficacy. "
        "Values < 1 μM are generally considered good candidates for drug development."
    ),
    "affinity": (
        "Affinity Probability is a predicted measure (0-1) of how likely the drug-protein complex "
        "exhibits strong binding. Higher values (closer to 1) suggest stronger predicted binding affinity "
        "between the drug molecule and the target protein. This metric is derived from machine learning "
        "models trained on experimental binding data."
    ),
    "conf": (
        "Confidence is a metric (0-1) from Boltz structure prediction indicating the overall quality "
        "and reliability of the predicted protein-ligand complex structure. Higher values indicate "
        "greater confidence in the structural prediction. Values > 0.7 are typically considered reliable."
    ),
    "ptm": (
        "pTM (predicted TM-score) is a metric (0-1) measuring the predicted accuracy of the overall "
        "protein fold compared to the true structure. Values > 0.5 indicate likely correct fold topology, "
        "while values > 0.8 suggest high-confidence accurate predictions. This metric comes from AlphaFold-based "
        "structure prediction methods."
    ),
    "iptm": (
        "ipTM (interface predicted TM-score) specifically measures the predicted accuracy of protein-ligand "
        "interface regions in the complex structure. Higher values (closer to 1) indicate greater confidence "
        "in the predicted binding site geometry and protein-ligand interactions. Values > 0.7 suggest "
        "high-quality interface predictions."
    ),
    "plddt": (
        "Average pLDDT (predicted Local Distance Difference Test) is a per-residue confidence score "
        "averaged across all residues in the structure. It ranges from 0-100, with values > 90 indicating "
        "very high confidence, 70-90 good confidence, 50-70 low confidence, and < 50 very low confidence. "
        "This metric helps identify which regions of the structure are reliable."
    ),
}
