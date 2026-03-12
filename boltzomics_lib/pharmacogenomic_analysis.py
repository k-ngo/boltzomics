"""
Pharmacogenomic Analysis Module for BoltzOmics

This module provides post-processing analysis of Boltz-2 predictions to extract
clinically actionable insights for pharmacogenomics applications.

Novel Contributions:
1. ΔpIC50 Analysis - Compute and classify mutation effects on drug binding
2. Confidence-Weighted Predictions - Use structural confidence for reliability
3. Ensemble Uncertainty Quantification - Leverage multi-model predictions
4. Mutation Effect Classification - Automatic clinical categorization
5. Binding Site Perturbation Analysis - Structural impact assessment

These analyses transform raw Boltz-2 outputs into clinically interpretable
pharmacogenomic insights.
"""

import numpy as np
import json
import os
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field, asdict
from enum import Enum
import math


class MutationEffect(Enum):
    """Classification of mutation effects on drug binding."""
    STRONG_SENSITIZING = "strong_sensitizing"      # ΔpIC50 > 1.0 (10x increase)
    MODERATE_SENSITIZING = "moderate_sensitizing"  # 0.5 < ΔpIC50 ≤ 1.0
    WEAK_SENSITIZING = "weak_sensitizing"          # 0.3 < ΔpIC50 ≤ 0.5
    NEUTRAL = "neutral"                            # -0.3 ≤ ΔpIC50 ≤ 0.3
    WEAK_RESISTANCE = "weak_resistance"            # -0.5 ≤ ΔpIC50 < -0.3
    MODERATE_RESISTANCE = "moderate_resistance"    # -1.0 ≤ ΔpIC50 < -0.5
    STRONG_RESISTANCE = "strong_resistance"        # ΔpIC50 < -1.0 (10x decrease)


class ConfidenceLevel(Enum):
    """Confidence levels for predictions."""
    HIGH = "high"           # ipTM > 0.8, pLDDT > 0.85
    MEDIUM = "medium"       # 0.6 < ipTM ≤ 0.8
    LOW = "low"             # ipTM ≤ 0.6
    UNRELIABLE = "unreliable"  # ipTM < 0.4


@dataclass
class PredictionUncertainty:
    """Uncertainty quantification for a prediction."""
    mean_pic50: float
    std_pic50: float
    confidence_interval_95: Tuple[float, float]
    ensemble_agreement: float  # 0-1, how well ensemble models agree
    structural_confidence: float  # Based on ipTM/pLDDT


@dataclass
class MutationAnalysis:
    """Complete analysis of a mutation's effect on drug binding."""
    mutation_name: str
    drug_name: str

    # Raw predictions
    wt_pic50: float
    mutant_pic50: float

    # Differential analysis
    delta_pic50: float
    fold_change: float  # 10^(ΔpIC50)

    # Classification
    effect_class: MutationEffect
    effect_description: str

    # Confidence metrics
    wt_confidence: float
    mutant_confidence: float
    prediction_reliability: ConfidenceLevel

    # Uncertainty
    uncertainty: Optional[PredictionUncertainty] = None

    # Clinical interpretation
    clinical_relevance: str = ""
    actionability_score: float = 0.0


@dataclass
class DrugSensitivityProfile:
    """Sensitivity profile of a protein variant across multiple drugs."""
    variant_name: str
    analyses: List[MutationAnalysis] = field(default_factory=list)

    # Aggregate metrics
    mean_delta_pic50: float = 0.0
    resistance_drugs: List[str] = field(default_factory=list)
    sensitizing_drugs: List[str] = field(default_factory=list)
    neutral_drugs: List[str] = field(default_factory=list)

    # Overall assessment
    overall_effect: str = ""
    drug_class_patterns: Dict[str, str] = field(default_factory=dict)


class PharmacogenomicAnalyzer:
    """
    Analyzes Boltz-2 predictions for pharmacogenomic insights.

    This class implements novel post-processing methods that transform
    raw structure predictions into clinically actionable information.
    """

    # Thresholds for effect classification (in pIC50 units)
    STRONG_THRESHOLD = 1.0      # 10-fold change
    MODERATE_THRESHOLD = 0.5    # ~3-fold change
    WEAK_THRESHOLD = 0.3        # ~2-fold change

    # Confidence thresholds
    HIGH_CONFIDENCE_IPTM = 0.8
    MEDIUM_CONFIDENCE_IPTM = 0.6
    LOW_CONFIDENCE_IPTM = 0.4

    def __init__(self):
        self.results_cache: Dict[str, Any] = {}

    def compute_delta_pic50(
        self,
        wt_affinity: Dict[str, float],
        mutant_affinity: Dict[str, float]
    ) -> Tuple[float, PredictionUncertainty]:
        """
        Compute ΔpIC50 with uncertainty quantification.

        Uses ensemble predictions from Boltz-2 to estimate uncertainty.
        Boltz-2 provides two ensemble models (value1, value2) which we
        use for uncertainty estimation.

        Args:
            wt_affinity: WT affinity dict from Boltz-2
            mutant_affinity: Mutant affinity dict from Boltz-2

        Returns:
            Tuple of (delta_pic50, uncertainty)
        """
        # Extract ensemble predictions
        # Boltz-2 outputs: affinity_pred_value (ensemble), value1, value2
        wt_values = [
            self._affinity_to_pic50(wt_affinity.get("affinity_pred_value", 0)),
            self._affinity_to_pic50(wt_affinity.get("affinity_pred_value1", 0)),
            self._affinity_to_pic50(wt_affinity.get("affinity_pred_value2", 0))
        ]

        mut_values = [
            self._affinity_to_pic50(mutant_affinity.get("affinity_pred_value", 0)),
            self._affinity_to_pic50(mutant_affinity.get("affinity_pred_value1", 0)),
            self._affinity_to_pic50(mutant_affinity.get("affinity_pred_value2", 0))
        ]

        # Compute delta for each ensemble member
        deltas = [m - w for w, m in zip(wt_values, mut_values)]

        mean_delta = np.mean(deltas)
        std_delta = np.std(deltas) if len(deltas) > 1 else 0.0

        # 95% confidence interval
        ci_95 = (
            mean_delta - 1.96 * std_delta,
            mean_delta + 1.96 * std_delta
        )

        # Ensemble agreement (inverse of coefficient of variation)
        ensemble_agreement = 1.0 - min(abs(std_delta / (abs(mean_delta) + 0.01)), 1.0)

        uncertainty = PredictionUncertainty(
            mean_pic50=mean_delta,
            std_pic50=std_delta,
            confidence_interval_95=ci_95,
            ensemble_agreement=ensemble_agreement,
            structural_confidence=0.0  # Set separately from confidence data
        )

        return mean_delta, uncertainty

    def _affinity_to_pic50(self, affinity_value: float) -> float:
        """
        Convert Boltz-2 affinity output to pIC50.

        Boltz-2 outputs log10(IC50) where IC50 is in μM.
        pIC50 = -log10(IC50 in M) = 6 - affinity_value
        """
        return 6.0 - affinity_value

    def classify_mutation_effect(self, delta_pic50: float) -> Tuple[MutationEffect, str]:
        """
        Classify the mutation effect based on ΔpIC50.

        Classification scheme:
        - Strong sensitizing: ΔpIC50 > 1.0 (>10-fold increase in affinity)
        - Moderate sensitizing: 0.5 < ΔpIC50 ≤ 1.0
        - Weak sensitizing: 0.3 < ΔpIC50 ≤ 0.5
        - Neutral: -0.3 ≤ ΔpIC50 ≤ 0.3
        - Weak resistance: -0.5 ≤ ΔpIC50 < -0.3
        - Moderate resistance: -1.0 ≤ ΔpIC50 < -0.5
        - Strong resistance: ΔpIC50 < -1.0 (>10-fold decrease)

        Returns:
            Tuple of (effect_class, description)
        """
        if delta_pic50 > self.STRONG_THRESHOLD:
            return MutationEffect.STRONG_SENSITIZING, \
                f"Strong sensitizing: {10**delta_pic50:.1f}-fold increase in binding affinity"
        elif delta_pic50 > self.MODERATE_THRESHOLD:
            return MutationEffect.MODERATE_SENSITIZING, \
                f"Moderate sensitizing: {10**delta_pic50:.1f}-fold increase in binding affinity"
        elif delta_pic50 > self.WEAK_THRESHOLD:
            return MutationEffect.WEAK_SENSITIZING, \
                f"Weak sensitizing: {10**delta_pic50:.1f}-fold increase in binding affinity"
        elif delta_pic50 >= -self.WEAK_THRESHOLD:
            return MutationEffect.NEUTRAL, \
                f"Neutral: No significant change in binding affinity ({10**abs(delta_pic50):.1f}-fold)"
        elif delta_pic50 >= -self.MODERATE_THRESHOLD:
            return MutationEffect.WEAK_RESISTANCE, \
                f"Weak resistance: {10**abs(delta_pic50):.1f}-fold decrease in binding affinity"
        elif delta_pic50 >= -self.STRONG_THRESHOLD:
            return MutationEffect.MODERATE_RESISTANCE, \
                f"Moderate resistance: {10**abs(delta_pic50):.1f}-fold decrease in binding affinity"
        else:
            return MutationEffect.STRONG_RESISTANCE, \
                f"Strong resistance: {10**abs(delta_pic50):.1f}-fold decrease in binding affinity"

    def assess_confidence(
        self,
        confidence_data: Dict[str, float]
    ) -> Tuple[ConfidenceLevel, float]:
        """
        Assess prediction confidence from Boltz-2 confidence metrics.

        Uses ligand_iptm (interface predicted TM score for ligand) as the
        primary indicator of binding prediction quality.

        Args:
            confidence_data: Confidence dict from Boltz-2

        Returns:
            Tuple of (confidence_level, numeric_score)
        """
        # Primary metric: ligand interface pTM
        ligand_iptm = confidence_data.get("ligand_iptm", 0.0)

        # Secondary metrics for composite score
        complex_plddt = confidence_data.get("complex_plddt", 0.0)
        iptm = confidence_data.get("iptm", 0.0)

        # Composite confidence score (weighted average)
        composite = 0.5 * ligand_iptm + 0.3 * iptm + 0.2 * complex_plddt

        if ligand_iptm > self.HIGH_CONFIDENCE_IPTM and complex_plddt > 0.85:
            return ConfidenceLevel.HIGH, composite
        elif ligand_iptm > self.MEDIUM_CONFIDENCE_IPTM:
            return ConfidenceLevel.MEDIUM, composite
        elif ligand_iptm > self.LOW_CONFIDENCE_IPTM:
            return ConfidenceLevel.LOW, composite
        else:
            return ConfidenceLevel.UNRELIABLE, composite

    def compute_clinical_actionability(
        self,
        effect: MutationEffect,
        confidence: ConfidenceLevel,
        delta_pic50: float,
        uncertainty: PredictionUncertainty
    ) -> Tuple[float, str]:
        """
        Compute clinical actionability score and interpretation.

        Actionability depends on:
        1. Magnitude of effect (|ΔpIC50|)
        2. Confidence in prediction
        3. Consistency of effect direction across ensemble

        Returns:
            Tuple of (actionability_score 0-1, clinical_interpretation)
        """
        # Base score from effect magnitude
        magnitude_score = min(abs(delta_pic50) / self.STRONG_THRESHOLD, 1.0)

        # Confidence modifier
        confidence_weights = {
            ConfidenceLevel.HIGH: 1.0,
            ConfidenceLevel.MEDIUM: 0.7,
            ConfidenceLevel.LOW: 0.4,
            ConfidenceLevel.UNRELIABLE: 0.1
        }
        confidence_modifier = confidence_weights.get(confidence, 0.5)

        # Ensemble agreement modifier
        agreement_modifier = uncertainty.ensemble_agreement if uncertainty else 0.5

        # Final actionability score
        actionability = magnitude_score * confidence_modifier * agreement_modifier

        # Clinical interpretation
        if actionability > 0.7:
            if effect in [MutationEffect.STRONG_RESISTANCE, MutationEffect.MODERATE_RESISTANCE]:
                interpretation = "HIGH PRIORITY: Strong evidence for drug resistance. Consider alternative therapy."
            elif effect in [MutationEffect.STRONG_SENSITIZING, MutationEffect.MODERATE_SENSITIZING]:
                interpretation = "HIGH PRIORITY: Strong evidence for increased drug sensitivity. Monitor for toxicity."
            else:
                interpretation = "Moderate clinical relevance. Standard monitoring recommended."
        elif actionability > 0.4:
            interpretation = "MODERATE PRIORITY: Consider in clinical decision-making with additional evidence."
        else:
            interpretation = "LOW PRIORITY: Insufficient evidence for clinical action. Further validation needed."

        return actionability, interpretation

    def analyze_mutation(
        self,
        mutation_name: str,
        drug_name: str,
        wt_affinity: Dict[str, float],
        mutant_affinity: Dict[str, float],
        wt_confidence: Dict[str, float],
        mutant_confidence: Dict[str, float]
    ) -> MutationAnalysis:
        """
        Perform complete analysis of a mutation's effect on drug binding.

        This is the main entry point for analyzing a single mutation-drug pair.

        Args:
            mutation_name: Name of the mutation (e.g., "K262R")
            drug_name: Name of the drug
            wt_affinity: WT affinity results from Boltz-2
            mutant_affinity: Mutant affinity results from Boltz-2
            wt_confidence: WT confidence metrics from Boltz-2
            mutant_confidence: Mutant confidence metrics from Boltz-2

        Returns:
            Complete MutationAnalysis object
        """
        # Compute ΔpIC50 with uncertainty
        delta_pic50, uncertainty = self.compute_delta_pic50(wt_affinity, mutant_affinity)

        # Get individual pIC50 values
        wt_pic50 = self._affinity_to_pic50(wt_affinity.get("affinity_pred_value", 0))
        mutant_pic50 = self._affinity_to_pic50(mutant_affinity.get("affinity_pred_value", 0))

        # Classify effect
        effect_class, effect_description = self.classify_mutation_effect(delta_pic50)

        # Assess confidence
        wt_conf_level, wt_conf_score = self.assess_confidence(wt_confidence)
        mut_conf_level, mut_conf_score = self.assess_confidence(mutant_confidence)

        # Use minimum confidence as overall reliability
        overall_conf_score = min(wt_conf_score, mut_conf_score)
        if overall_conf_score > 0.75:
            reliability = ConfidenceLevel.HIGH
        elif overall_conf_score > 0.55:
            reliability = ConfidenceLevel.MEDIUM
        elif overall_conf_score > 0.35:
            reliability = ConfidenceLevel.LOW
        else:
            reliability = ConfidenceLevel.UNRELIABLE

        # Update uncertainty with structural confidence
        uncertainty.structural_confidence = overall_conf_score

        # Compute clinical actionability
        actionability, clinical_relevance = self.compute_clinical_actionability(
            effect_class, reliability, delta_pic50, uncertainty
        )

        return MutationAnalysis(
            mutation_name=mutation_name,
            drug_name=drug_name,
            wt_pic50=wt_pic50,
            mutant_pic50=mutant_pic50,
            delta_pic50=delta_pic50,
            fold_change=10 ** delta_pic50,
            effect_class=effect_class,
            effect_description=effect_description,
            wt_confidence=wt_conf_score,
            mutant_confidence=mut_conf_score,
            prediction_reliability=reliability,
            uncertainty=uncertainty,
            clinical_relevance=clinical_relevance,
            actionability_score=actionability
        )

    def create_sensitivity_profile(
        self,
        variant_name: str,
        analyses: List[MutationAnalysis]
    ) -> DrugSensitivityProfile:
        """
        Create a drug sensitivity profile for a protein variant.

        Aggregates analyses across multiple drugs to identify patterns.
        """
        profile = DrugSensitivityProfile(
            variant_name=variant_name,
            analyses=analyses
        )

        if not analyses:
            return profile

        # Aggregate metrics
        deltas = [a.delta_pic50 for a in analyses]
        profile.mean_delta_pic50 = np.mean(deltas)

        # Categorize drugs by effect
        for analysis in analyses:
            if analysis.effect_class in [MutationEffect.STRONG_RESISTANCE,
                                         MutationEffect.MODERATE_RESISTANCE]:
                profile.resistance_drugs.append(analysis.drug_name)
            elif analysis.effect_class in [MutationEffect.STRONG_SENSITIZING,
                                           MutationEffect.MODERATE_SENSITIZING]:
                profile.sensitizing_drugs.append(analysis.drug_name)
            else:
                profile.neutral_drugs.append(analysis.drug_name)

        # Overall assessment
        n_resist = len(profile.resistance_drugs)
        n_sensit = len(profile.sensitizing_drugs)
        n_total = len(analyses)

        if n_resist > n_total * 0.5:
            profile.overall_effect = f"Multi-drug resistance variant ({n_resist}/{n_total} drugs affected)"
        elif n_sensit > n_total * 0.5:
            profile.overall_effect = f"Multi-drug sensitizing variant ({n_sensit}/{n_total} drugs affected)"
        elif n_resist > 0 or n_sensit > 0:
            profile.overall_effect = f"Drug-specific effects (R:{n_resist}, S:{n_sensit}, N:{len(profile.neutral_drugs)})"
        else:
            profile.overall_effect = "Neutral variant - no significant drug interaction changes"

        return profile

    def generate_report(
        self,
        analyses: List[MutationAnalysis],
        output_format: str = "dict"
    ) -> Any:
        """
        Generate a comprehensive pharmacogenomic report.

        Args:
            analyses: List of mutation analyses
            output_format: "dict", "json", or "markdown"

        Returns:
            Report in specified format
        """
        # Group by variant
        variants: Dict[str, List[MutationAnalysis]] = {}
        for analysis in analyses:
            if analysis.mutation_name not in variants:
                variants[analysis.mutation_name] = []
            variants[analysis.mutation_name].append(analysis)

        # Create profiles
        profiles = [
            self.create_sensitivity_profile(name, analyses)
            for name, analyses in variants.items()
        ]

        # Summary statistics
        all_deltas = [a.delta_pic50 for a in analyses]
        high_impact = [a for a in analyses if abs(a.delta_pic50) > self.MODERATE_THRESHOLD]
        high_confidence = [a for a in analyses if a.prediction_reliability == ConfidenceLevel.HIGH]

        report = {
            "summary": {
                "total_predictions": len(analyses),
                "variants_analyzed": len(variants),
                "drugs_tested": len(set(a.drug_name for a in analyses)),
                "high_impact_interactions": len(high_impact),
                "high_confidence_predictions": len(high_confidence),
                "mean_delta_pic50": float(np.mean(all_deltas)) if all_deltas else 0.0,
                "std_delta_pic50": float(np.std(all_deltas)) if all_deltas else 0.0
            },
            "variant_profiles": [
                {
                    "variant": p.variant_name,
                    "overall_effect": p.overall_effect,
                    "resistance_drugs": p.resistance_drugs,
                    "sensitizing_drugs": p.sensitizing_drugs,
                    "neutral_drugs": p.neutral_drugs,
                    "mean_delta_pic50": p.mean_delta_pic50
                }
                for p in profiles
            ],
            "high_priority_findings": [
                {
                    "mutation": a.mutation_name,
                    "drug": a.drug_name,
                    "effect": a.effect_class.value,
                    "delta_pic50": a.delta_pic50,
                    "fold_change": a.fold_change,
                    "confidence": a.prediction_reliability.value,
                    "clinical_relevance": a.clinical_relevance,
                    "actionability": a.actionability_score
                }
                for a in sorted(high_impact, key=lambda x: -abs(x.delta_pic50))
            ],
            "all_analyses": [
                {
                    "mutation": a.mutation_name,
                    "drug": a.drug_name,
                    "wt_pic50": a.wt_pic50,
                    "mutant_pic50": a.mutant_pic50,
                    "delta_pic50": a.delta_pic50,
                    "fold_change": a.fold_change,
                    "effect": a.effect_class.value,
                    "confidence": a.prediction_reliability.value,
                    "actionability": a.actionability_score
                }
                for a in analyses
            ]
        }

        if output_format == "json":
            return json.dumps(report, indent=2)
        elif output_format == "markdown":
            return self._report_to_markdown(report)
        else:
            return report

    def _report_to_markdown(self, report: Dict) -> str:
        """Convert report to markdown format."""
        lines = [
            "# Pharmacogenomic Analysis Report",
            "",
            "## Summary",
            f"- **Total predictions:** {report['summary']['total_predictions']}",
            f"- **Variants analyzed:** {report['summary']['variants_analyzed']}",
            f"- **Drugs tested:** {report['summary']['drugs_tested']}",
            f"- **High-impact interactions:** {report['summary']['high_impact_interactions']}",
            f"- **High-confidence predictions:** {report['summary']['high_confidence_predictions']}",
            f"- **Mean ΔpIC50:** {report['summary']['mean_delta_pic50']:.3f} ± {report['summary']['std_delta_pic50']:.3f}",
            "",
            "## Variant Profiles",
            ""
        ]

        for profile in report["variant_profiles"]:
            lines.extend([
                f"### {profile['variant']}",
                f"**Overall effect:** {profile['overall_effect']}",
                f"- Resistance drugs: {', '.join(profile['resistance_drugs']) or 'None'}",
                f"- Sensitizing drugs: {', '.join(profile['sensitizing_drugs']) or 'None'}",
                f"- Neutral drugs: {', '.join(profile['neutral_drugs']) or 'None'}",
                ""
            ])

        if report["high_priority_findings"]:
            lines.extend([
                "## High-Priority Findings",
                "",
                "| Mutation | Drug | Effect | ΔpIC50 | Fold Change | Confidence | Actionability |",
                "|----------|------|--------|--------|-------------|------------|---------------|"
            ])

            for finding in report["high_priority_findings"][:10]:  # Top 10
                lines.append(
                    f"| {finding['mutation']} | {finding['drug']} | {finding['effect']} | "
                    f"{finding['delta_pic50']:.2f} | {finding['fold_change']:.1f}x | "
                    f"{finding['confidence']} | {finding['actionability']:.2f} |"
                )

        return "\n".join(lines)


# Convenience functions for integration
def analyze_screening_results(
    results: List[Dict[str, Any]],
    wt_name: str = "WT"
) -> Dict[str, Any]:
    """
    Analyze a list of screening results from BoltzOmics.

    Args:
        results: List of result dicts from screening
        wt_name: Name of the wild-type protein

    Returns:
        Pharmacogenomic analysis report
    """
    analyzer = PharmacogenomicAnalyzer()

    # Group results by drug to pair WT with mutants
    by_drug: Dict[str, Dict[str, Dict]] = {}
    for result in results:
        drug = result.get("drug_name", result.get("ligand_name", "Unknown"))
        protein = result.get("protein_name", "Unknown")

        if drug not in by_drug:
            by_drug[drug] = {}
        by_drug[drug][protein] = result

    # Analyze each mutation-drug pair
    analyses = []
    for drug, proteins in by_drug.items():
        if wt_name not in proteins:
            continue

        wt_result = proteins[wt_name]
        wt_affinity = {
            "affinity_pred_value": wt_result.get("affinity_pred_value", 0),
            "affinity_pred_value1": wt_result.get("affinity_pred_value1", 0),
            "affinity_pred_value2": wt_result.get("affinity_pred_value2", 0)
        }
        wt_confidence = {
            "ligand_iptm": wt_result.get("ligand_iptm", wt_result.get("iptm", 0)),
            "iptm": wt_result.get("iptm", 0),
            "complex_plddt": wt_result.get("complex_plddt", wt_result.get("avg_plddt", 0) / 100)
        }

        for protein_name, result in proteins.items():
            if protein_name == wt_name:
                continue

            mut_affinity = {
                "affinity_pred_value": result.get("affinity_pred_value", 0),
                "affinity_pred_value1": result.get("affinity_pred_value1", 0),
                "affinity_pred_value2": result.get("affinity_pred_value2", 0)
            }
            mut_confidence = {
                "ligand_iptm": result.get("ligand_iptm", result.get("iptm", 0)),
                "iptm": result.get("iptm", 0),
                "complex_plddt": result.get("complex_plddt", result.get("avg_plddt", 0) / 100)
            }

            analysis = analyzer.analyze_mutation(
                mutation_name=protein_name,
                drug_name=drug,
                wt_affinity=wt_affinity,
                mutant_affinity=mut_affinity,
                wt_confidence=wt_confidence,
                mutant_confidence=mut_confidence
            )
            analyses.append(analysis)

    return analyzer.generate_report(analyses)


def quick_delta_analysis(
    wt_pic50: float,
    mutant_pic50: float,
    wt_confidence: float = 0.8,
    mutant_confidence: float = 0.8
) -> Dict[str, Any]:
    """
    Quick ΔpIC50 analysis for a single prediction pair.

    Args:
        wt_pic50: Wild-type pIC50
        mutant_pic50: Mutant pIC50
        wt_confidence: WT prediction confidence (0-1)
        mutant_confidence: Mutant prediction confidence (0-1)

    Returns:
        Analysis summary dict
    """
    analyzer = PharmacogenomicAnalyzer()

    delta = mutant_pic50 - wt_pic50
    effect_class, description = analyzer.classify_mutation_effect(delta)

    return {
        "delta_pic50": delta,
        "fold_change": 10 ** delta,
        "effect_class": effect_class.value,
        "description": description,
        "confidence": min(wt_confidence, mutant_confidence)
    }
