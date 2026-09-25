#!/usr/bin/env python3
"""Code-first submission and management for BoltzOmics screening jobs.

Job inputs live in a YAML or JSON manifest. Protein FASTA and ligand CSV files
can be referenced from that manifest, keeping large screening panels out of
the command itself.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml


CLI_STATE_DIRECTORY = "_cli_job_state"

# These are the public per-job settings consumed by BoltzOmics' screening
# executor. Nested values use the same structures as the Streamlit controls.
SUPPORTED_SETTINGS = {
    "use_gpu",
    "accelerator",
    "devices",
    "cuda_visible_devices",
    "queue_gpu_devices",
    "preprocessing_threads",
    "override",
    "recycling_steps",
    "sampling_steps",
    "diffusion_samples",
    "max_parallel_samples",
    "step_scale",
    "affinity_mw_correction",
    "affinity_consensus_enabled",
    "affinity_consensus_mode",
    "affinity_consensus_weight_floor",
    "affinity_consensus_entropy_alpha",
    "external_boltz_patch_enabled",
    "external_boltz_patch_mode",
    "external_boltz_patch_weight_floor",
    "external_boltz_patch_entropy_alpha",
    "external_boltz_patch_uncertainty_penalty",
    "external_boltz_patch_min_confidence",
    "affinity_multisampling_enabled",
    "affinity_multisampling_profiles",
    "affinity_multisampling_settings",
    "affinity_multisampling_refinement_steps",
    "affinity_multisampling_aggregate_mode",
    "affinity_multisampling_apply_aggregate",
    "affinity_multisampling_early_stop_enabled",
    "affinity_multisampling_early_stop_min_points",
    "affinity_multisampling_early_stop_delta",
    "affinity_multisampling_early_stop_std",
    "affinity_multisampling_early_stop_patience",
    "affinity_multisampling_robust_outlier_filter",
    "affinity_multisampling_robust_outlier_zmax",
    "affinity_multisampling_bootstrap_samples",
    "confidence_target",
    "max_msa_seqs",
    "sampling_steps_affinity",
    "diffusion_samples_affinity",
    "subsample_msa",
    "num_subsampled_msa",
    "use_potentials",
    "method",
    "mutation_steering_config",
    "template_cif_path",
    "template_options",
    "boltz_runtime_options",
    "binding_pocket_constraints",
    "distance_constraints",
    "cofactor_info",
    "ptm_modifications",
    "prediction_timeout_seconds",
    "enable_retries",
    "max_retry_attempts",
    "retry_delay_base",
    "enable_msa_cache",
    "prediction_replicas",
}


class JobSpecError(ValueError):
    """Raised when a job manifest or one of its input files is invalid."""


def _application():
    """Import the existing BoltzOmics prediction implementation on demand."""
    import boltzomics

    return boltzomics


def _manifest_file_path(raw_path: Any, manifest_dir: Path, label: str) -> Path:
    if not isinstance(raw_path, (str, os.PathLike)) or not str(raw_path).strip():
        raise JobSpecError(f"{label} must be a non-empty file path.")
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = manifest_dir / path
    path = path.resolve()
    if not path.is_file():
        raise JobSpecError(f"{label} file does not exist: {path}")
    return path


def _validate_project_name(name: Any) -> str:
    project_name = str(name or "").strip()
    if not project_name:
        raise JobSpecError("project_name is required.")
    if project_name in {".", ".."} or "/" in project_name or "\\" in project_name:
        raise JobSpecError("project_name must be a plain folder name without path separators.")
    return project_name


def _manifest_bool(value: Any, label: str) -> bool:
    if not isinstance(value, bool):
        raise JobSpecError(f"{label} must be true or false.")
    return value


def _load_proteins(raw: Any, manifest_dir: Path, app) -> List[Tuple[str, str]]:
    if isinstance(raw, dict) and ("fasta_file" in raw or "fasta" in raw):
        fasta_path = _manifest_file_path(raw.get("fasta_file", raw.get("fasta")), manifest_dir, "proteins.fasta")
        records = app.parse_fasta_sequences(fasta_path.read_text(encoding="utf-8"))
    elif isinstance(raw, str):
        fasta_path = _manifest_file_path(raw, manifest_dir, "proteins FASTA")
        records = app.parse_fasta_sequences(fasta_path.read_text(encoding="utf-8"))
    elif isinstance(raw, list):
        records = []
        for index, item in enumerate(raw, 1):
            if not isinstance(item, dict):
                raise JobSpecError(f"proteins[{index - 1}] must contain name and sequence fields.")
            name = str(item.get("name") or f"Protein_{index}").strip()
            sequence = str(item.get("sequence") or "").strip()
            records.append((name, sequence))
    else:
        raise JobSpecError("proteins must be a FASTA file reference or a list of {name, sequence} records.")

    if not records:
        raise JobSpecError("No protein sequences were found.")

    proteins: List[Tuple[str, str]] = []
    used_names = set()
    for index, (raw_name, raw_sequence) in enumerate(records, 1):
        name = str(raw_name or f"Protein_{index}").strip()
        if not name:
            raise JobSpecError(f"Protein {index} has an empty name.")
        if name.casefold() in used_names:
            raise JobSpecError(f"Protein names must be unique; repeated name: {name!r}.")
        used_names.add(name.casefold())
        valid, error, _chains, sequence = app.validate_protein_sequence(str(raw_sequence or ""))
        if not valid:
            raise JobSpecError(f"Invalid sequence for protein {name!r}: {error}")
        proteins.append((name, sequence))
    return proteins


def _load_ligands(raw: Any, manifest_dir: Path, app) -> List[Tuple[str, str]]:
    if raw is None:
        return []
    if isinstance(raw, dict) and ("csv_file" in raw or "csv" in raw):
        csv_path = _manifest_file_path(raw.get("csv_file", raw.get("csv")), manifest_dir, "ligands.csv")
        name_column = str(raw.get("name_column", "name"))
        smiles_column = str(raw.get("smiles_column", "smiles"))
        delimiter = "\t" if csv_path.suffix.lower() == ".tsv" else ","
        with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            if not reader.fieldnames:
                raise JobSpecError(f"Ligand file has no header row: {csv_path}")
            columns = {str(column).strip().casefold(): column for column in reader.fieldnames if column}
            actual_name_column = columns.get(name_column.casefold())
            actual_smiles_column = columns.get(smiles_column.casefold())
            if actual_smiles_column is None:
                raise JobSpecError(
                    f"Ligand file {csv_path} is missing the configured SMILES column {smiles_column!r}."
                )
            records = []
            for index, row in enumerate(reader, 1):
                name = row.get(actual_name_column, "") if actual_name_column else f"Drug_{index}"
                smiles = row.get(actual_smiles_column, "")
                if not str(name or "").strip() and not str(smiles or "").strip():
                    continue
                records.append((str(name or f"Drug_{index}").strip(), str(smiles or "").strip()))
    elif isinstance(raw, str):
        csv_path = _manifest_file_path(raw, manifest_dir, "ligands CSV")
        return _load_ligands({"csv_file": str(csv_path)}, manifest_dir, app)
    elif isinstance(raw, list):
        records = []
        for index, item in enumerate(raw, 1):
            if not isinstance(item, dict):
                raise JobSpecError(f"ligands[{index - 1}] must contain name and smiles fields.")
            records.append((str(item.get("name") or f"Drug_{index}").strip(), str(item.get("smiles") or "").strip()))
    else:
        raise JobSpecError("ligands must be a CSV file reference or a list of {name, smiles} records.")

    ligands: List[Tuple[str, str]] = []
    used_names = set()
    for index, (raw_name, smiles) in enumerate(records, 1):
        name = str(raw_name or f"Drug_{index}").strip()
        if not name:
            raise JobSpecError(f"Ligand {index} has an empty name.")
        if name.casefold() in used_names:
            raise JobSpecError(f"Ligand names must be unique; repeated name: {name!r}.")
        used_names.add(name.casefold())
        if not smiles:
            raise JobSpecError(f"Ligand {name!r} has an empty SMILES string.")
        if not app.validate_smiles(smiles):
            raise JobSpecError(f"Invalid SMILES for ligand {name!r}: {smiles!r}")
        ligands.append((name, smiles))
    return ligands


def _load_filter(raw: Any, proteins: Sequence[Tuple[str, str]], ligands: Sequence[Tuple[str, str]]) -> Optional[Dict[str, Any]]:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise JobSpecError("protein_drug_filter must be an object with enabled and pairs fields.")
    pairs = raw.get("pairs", [])
    if not isinstance(pairs, list):
        raise JobSpecError("protein_drug_filter.pairs must be a list.")
    normalized_pairs = []
    protein_names = {name for name, _ in proteins}
    ligand_names = {name for name, _ in ligands}
    for index, pair in enumerate(pairs, 1):
        if isinstance(pair, dict):
            protein_name = str(pair.get("protein") or "").strip()
            ligand_name = str(pair.get("ligand") or "").strip()
        elif isinstance(pair, (list, tuple)) and len(pair) == 2:
            protein_name, ligand_name = (str(value).strip() for value in pair)
        else:
            raise JobSpecError(f"protein_drug_filter.pairs[{index - 1}] must be [protein, ligand].")
        if protein_name not in protein_names:
            raise JobSpecError(f"Filter references unknown protein {protein_name!r}.")
        if ligand_name not in ligand_names:
            raise JobSpecError(f"Filter references unknown ligand {ligand_name!r}.")
        normalized_pairs.append((protein_name, ligand_name))
    enabled = _manifest_bool(raw.get("enabled", True), "protein_drug_filter.enabled")
    return {"enabled": enabled, "pairs": normalized_pairs}


def load_job_spec(spec_path: os.PathLike[str] | str) -> Dict[str, Any]:
    """Read and validate a YAML/JSON manifest and the files it references."""
    manifest_path = Path(spec_path).expanduser().resolve()
    if not manifest_path.is_file():
        raise JobSpecError(f"Job manifest does not exist: {manifest_path}")
    try:
        with manifest_path.open("r", encoding="utf-8") as handle:
            raw_spec = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as exc:
        raise JobSpecError(f"Could not read job manifest {manifest_path}: {exc}") from exc
    if not isinstance(raw_spec, dict):
        raise JobSpecError("Job manifest must be a YAML or JSON object.")

    app = _application()
    project_name = _validate_project_name(raw_spec.get("project_name"))
    structure_only = _manifest_bool(raw_spec.get("structure_only", False), "structure_only")
    proteins = _load_proteins(raw_spec.get("proteins"), manifest_path.parent, app)
    ligands = _load_ligands(raw_spec.get("ligands"), manifest_path.parent, app)
    if not structure_only and not ligands:
        raise JobSpecError("Provide at least one ligand, or set structure_only: true.")

    gpu_mode = str(raw_spec.get("gpu_mode", "auto")).strip().lower()
    if gpu_mode not in {"auto", "single"}:
        raise JobSpecError("gpu_mode must be either 'auto' or 'single'.")

    raw_settings = raw_spec.get("settings", raw_spec.get("parameters", {}))
    if not isinstance(raw_settings, dict):
        raise JobSpecError("settings must be an object of BoltzOmics prediction options.")
    unknown_settings = sorted(set(raw_settings) - SUPPORTED_SETTINGS)
    if unknown_settings:
        raise JobSpecError("Unsupported setting(s): " + ", ".join(unknown_settings))
    settings = dict(raw_settings)
    mapping_settings = (
        "mutation_steering_config",
        "template_options",
        "boltz_runtime_options",
        "binding_pocket_constraints",
    )
    for key in mapping_settings:
        if key in settings and settings[key] is not None and not isinstance(settings[key], dict):
            raise JobSpecError(f"settings.{key} must be an object.")
    list_settings = (
        "queue_gpu_devices",
        "affinity_multisampling_profiles",
        "affinity_multisampling_settings",
        "affinity_multisampling_refinement_steps",
        "distance_constraints",
    )
    for key in list_settings:
        if key in settings and settings[key] is not None and not isinstance(settings[key], list):
            raise JobSpecError(f"settings.{key} must be a list.")
    if "queue_gpu_devices" in settings and settings["queue_gpu_devices"] is not None:
        device_ids = settings["queue_gpu_devices"]
        if not device_ids or any(isinstance(device, bool) or not str(device).strip() for device in device_ids):
            raise JobSpecError("settings.queue_gpu_devices must contain non-empty GPU device IDs.")
        device_ids = [str(device).strip() for device in device_ids]
        if len(set(device_ids)) != len(device_ids):
            raise JobSpecError("settings.queue_gpu_devices cannot contain duplicate device IDs.")
        settings["queue_gpu_devices"] = device_ids
        if gpu_mode == "single":
            raise JobSpecError("Do not set queue_gpu_devices when gpu_mode is 'single'.")
    if "cofactor_info" in settings and settings["cofactor_info"] is not None:
        cofactor_info = settings["cofactor_info"]
        if not isinstance(cofactor_info, (dict, list)):
            raise JobSpecError("settings.cofactor_info must be an object or a list of objects.")
        if isinstance(cofactor_info, list) and any(not isinstance(item, dict) for item in cofactor_info):
            raise JobSpecError("Each settings.cofactor_info entry must be an object.")
    if "ptm_modifications" in settings and settings["ptm_modifications"] is not None:
        if not isinstance(settings["ptm_modifications"], dict):
            raise JobSpecError("settings.ptm_modifications must be an object.")
    profiles = settings.get("affinity_multisampling_profiles")
    if profiles is not None and any(not isinstance(profile, dict) for profile in profiles):
        raise JobSpecError("Each affinity_multisampling_profiles entry must be an object.")
    if profiles is not None:
        for index, profile in enumerate(profiles, 1):
            if not {"sampling_steps_affinity", "diffusion_samples_affinity"}.issubset(profile):
                raise JobSpecError(
                    f"affinity_multisampling_profiles entry {index} needs sampling_steps_affinity "
                    "and diffusion_samples_affinity."
                )
    if settings.get("distance_constraints") is not None:
        try:
            settings["distance_constraints"] = app.normalize_distance_constraints(
                settings["distance_constraints"]
            )
        except ValueError as exc:
            raise JobSpecError(f"settings.distance_constraints is invalid: {exc}") from exc
    if "template_cif_path" in settings and settings["template_cif_path"]:
        settings["template_cif_path"] = str(
            _manifest_file_path(settings["template_cif_path"], manifest_path.parent, "settings.template_cif_path")
        )
    pocket_forced = bool((settings.get("binding_pocket_constraints") or {}).get("force", False))
    distance_forced = any(
        constraint.get("force", False) for constraint in settings.get("distance_constraints", [])
    )
    template_forced = bool((settings.get("template_options") or {}).get("force", False))
    steering = settings.get("mutation_steering_config") or {}
    steering_forced = bool(steering.get("enabled", False) and steering.get("use_potentials", False))
    settings["use_potentials"] = bool(
        settings.get("use_potentials", False)
        or pocket_forced
        or distance_forced
        or template_forced
        or steering_forced
    )

    default_accelerator = "cpu" if settings.get("use_gpu") is False else "gpu"
    accelerator = str(settings.get("accelerator", default_accelerator)).strip().lower()
    if accelerator not in {"gpu", "cpu"}:
        raise JobSpecError("settings.accelerator must be 'gpu' or 'cpu'.")
    use_gpu = settings.get("use_gpu", accelerator == "gpu")
    if not isinstance(use_gpu, bool):
        raise JobSpecError("settings.use_gpu must be true or false.")
    if use_gpu != (accelerator == "gpu"):
        raise JobSpecError("settings.use_gpu and settings.accelerator must select the same device type.")
    settings["accelerator"] = accelerator
    settings["use_gpu"] = use_gpu
    if settings.get("queue_gpu_devices") and not use_gpu:
        raise JobSpecError("settings.queue_gpu_devices requires GPU execution.")

    # Auto mode distributes queue workers across every detected CUDA GPU. An
    # explicit queue_gpu_devices list takes precedence over discovery.
    if (
        gpu_mode == "auto"
        and use_gpu
        and "queue_gpu_devices" not in settings
    ):
        discovered_gpus = app.discover_gpu_devices()
        if len(discovered_gpus) > 1:
            settings["queue_gpu_devices"] = [str(device_id) for device_id, _name in discovered_gpus]

    use_existing_results = _manifest_bool(raw_spec.get("use_existing_results", True), "use_existing_results")
    settings.setdefault("override", not use_existing_results)
    try:
        raw_workers = raw_spec.get("workers")
        suggested_workers = len(settings.get("queue_gpu_devices") or []) or 1
        workers = int(suggested_workers if raw_workers is None else raw_workers)
    except (TypeError, ValueError) as exc:
        raise JobSpecError("workers must be a positive integer.") from exc
    if workers < 1:
        raise JobSpecError("workers must be a positive integer.")
    if raw_workers is not None and settings.get("queue_gpu_devices"):
        settings["queue_gpu_devices"] = settings["queue_gpu_devices"][:workers]

    return {
        "manifest_path": str(manifest_path),
        "project_name": project_name,
        "structure_only": structure_only,
        "proteins": proteins,
        "ligands": ligands,
        "use_existing_results": use_existing_results,
        "protein_drug_filter": _load_filter(raw_spec.get("protein_drug_filter"), proteins, ligands),
        "settings": settings,
        "workers": workers,
        "gpu_mode": gpu_mode,
    }


def _new_manager(app):
    # Keep command-line worker state separate from the Streamlit queue while
    # writing predictions and project results into the same results directory.
    state_dir = os.path.join(app.RESULTS_DIR, CLI_STATE_DIRECTORY)
    return app.ScreeningJobManager(state_dir)


def _queue_counts(jobs: Iterable[Any]) -> Dict[str, int]:
    counts = {"pending": 0, "running": 0, "success": 0, "failed": 0, "cancelled": 0, "total": 0}
    for job in jobs:
        status = str(getattr(job, "status", ""))
        counts["total"] += 1
        if status in counts:
            counts[status] += 1
    return counts


def _suggest_worker_count(jobs: Iterable[Any], minimum: int = 1) -> int:
    counts = [max(1, int(minimum))]
    for job in jobs:
        if getattr(job, "status", "pending") not in {"pending", "running"}:
            continue
        params = getattr(job, "parameters", {}) or {}
        if params.get("use_gpu", True) and str(params.get("accelerator", "gpu")).lower() == "gpu":
            devices = params.get("queue_gpu_devices")
            if isinstance(devices, list) and devices:
                counts.append(len(devices))
    return max(counts)


def _queue_manifest(spec: Dict[str, Any], app, manager) -> Dict[str, Any]:
    jobs, cached_results, preparation = app.prepare_screening_jobs(
        protein_sequences=spec["proteins"],
        drug_smiles=spec["ligands"],
        project_name=spec["project_name"],
        structure_only=spec["structure_only"],
        use_existing_results=spec["use_existing_results"],
        protein_drug_filter=spec["protein_drug_filter"],
        shared_params=spec["settings"],
        manager=manager,
    )
    enqueued = manager.enqueue_jobs(jobs)
    app.ensure_project_metadata(spec["project_name"], app.RESULTS_DIR)
    return {
        "project_name": spec["project_name"],
        "manifest": spec["manifest_path"],
        "jobs_prepared": len(jobs),
        "jobs_enqueued": len(enqueued),
        "job_ids": [job.job_id for job in enqueued],
        "cached_results": len(cached_results),
        "skipped": preparation["skipped"],
        "duplicates": preparation["duplicate_jobs"],
        "warnings": preparation["warnings"],
        "queue": _queue_counts(manager.get_project_jobs(spec["project_name"])),
    }


def submit_spec(spec_path: os.PathLike[str] | str) -> Dict[str, Any]:
    """Queue jobs from a manifest without starting prediction workers."""
    spec = load_job_spec(spec_path)
    app = _application()
    manager = _new_manager(app)
    report = _queue_manifest(spec, app, manager)
    report["workers_started"] = False
    report["workers"] = spec["workers"]
    report["gpu_mode"] = spec["gpu_mode"]
    return report


def _start_workers(manager, app, worker_count: int) -> None:
    manager.set_worker_count(worker_count)
    manager.register_executor(app.execute_screening_job, failure_handler=app._persist_job_failure)


def run_spec(spec_path: os.PathLike[str] | str, workers: Optional[int] = None) -> Dict[str, Any]:
    """Queue and run a manifest, returning after the CLI queue is idle."""
    spec = load_job_spec(spec_path)
    app = _application()
    worker_count = int(spec["workers"] if workers is None else workers)
    if worker_count < 1:
        raise JobSpecError("workers must be a positive integer.")
    manager = _new_manager(app)
    report = _queue_manifest(spec, app, manager)
    if workers is None:
        with manager.lock:
            worker_count = _suggest_worker_count(manager.jobs.values(), minimum=spec["workers"])
    _start_workers(manager, app, worker_count)
    try:
        while True:
            project_counts = _queue_counts(manager.get_project_jobs(spec["project_name"]))
            with manager.lock:
                all_jobs = list(manager.jobs.values())
            all_counts = _queue_counts(all_jobs)
            if all_counts["pending"] == 0 and all_counts["running"] == 0:
                report["queue"] = project_counts
                report["workers_started"] = True
                report["workers"] = worker_count
                report["gpu_mode"] = spec["gpu_mode"]
                report["exit_code"] = 1 if project_counts["failed"] else 0
                return report
            time.sleep(0.5)
    finally:
        manager.shutdown(purge_pending=False)


def worker_loop(workers: Optional[int] = None, poll_seconds: float = 1.0) -> Dict[str, int]:
    """Process the persisted CLI queue and exit after it drains."""
    if workers is not None and workers < 1:
        raise JobSpecError("workers must be a positive integer.")
    if poll_seconds <= 0:
        raise JobSpecError("poll_seconds must be greater than zero.")
    app = _application()
    manager = _new_manager(app)
    with manager.lock:
        worker_count = _suggest_worker_count(manager.jobs.values()) if workers is None else workers
    _start_workers(manager, app, worker_count)
    try:
        while True:
            with manager.lock:
                counts = _queue_counts(list(manager.jobs.values()))
            if counts["pending"] == 0 and counts["running"] == 0:
                return counts
            time.sleep(poll_seconds)
    finally:
        manager.shutdown(purge_pending=False)


def _read_cli_state(app) -> List[Dict[str, Any]]:
    state_path = os.path.join(app.RESULTS_DIR, CLI_STATE_DIRECTORY, "job_state.json")
    if not os.path.isfile(state_path):
        return []
    try:
        with open(state_path, "r", encoding="utf-8") as handle:
            state = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise JobSpecError(f"Could not read CLI queue state: {exc}") from exc
    jobs = state.get("jobs", []) if isinstance(state, dict) else []
    return [job for job in jobs if isinstance(job, dict)]


def status_report(project_name: Optional[str] = None) -> Dict[str, Any]:
    """Read persisted CLI queue state without changing it."""
    app = _application()
    jobs = _read_cli_state(app)
    if project_name:
        jobs = [job for job in jobs if job.get("project_name") == project_name]
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for job in jobs:
        grouped.setdefault(str(job.get("project_name", "")), []).append(job)

    def summarize(project_jobs: List[Dict[str, Any]]) -> Dict[str, Any]:
        counts = {"pending": 0, "running": 0, "success": 0, "failed": 0, "cancelled": 0, "total": len(project_jobs)}
        for job in project_jobs:
            status = str(job.get("status", ""))
            if status in counts and status != "total":
                counts[status] += 1
        counts["jobs"] = [
            {
                "job_id": job.get("job_id"),
                "protein": job.get("protein_name"),
                "ligand": job.get("drug_name") or None,
                "status": job.get("status"),
                "retries": job.get("retries", 0),
                "error": job.get("error"),
            }
            for job in project_jobs
        ]
        return counts

    if project_name:
        return {"project_name": project_name, **summarize(grouped.get(project_name, []))}
    return {"projects": {name: summarize(items) for name, items in sorted(grouped.items())}}


def wait_for_project(project_name: str, timeout: Optional[float] = None, poll_seconds: float = 2.0) -> Dict[str, Any]:
    """Wait for a CLI worker process to finish a project's pending jobs."""
    if poll_seconds <= 0:
        raise JobSpecError("poll_seconds must be greater than zero.")
    if timeout is not None and timeout < 0:
        raise JobSpecError("timeout cannot be negative.")
    start = time.monotonic()
    while True:
        report = status_report(project_name)
        if report["pending"] == 0 and report["running"] == 0:
            return report
        if timeout is not None and time.monotonic() - start >= timeout:
            raise TimeoutError(f"Timed out waiting for project {project_name!r}.")
        time.sleep(poll_seconds)


def retry_project(project_name: str) -> Dict[str, Any]:
    app = _application()
    manager = _new_manager(app)
    requeued = manager.retry_failed_jobs(project_name)
    return {"project_name": project_name, "requeued": requeued, "queue": _queue_counts(manager.get_project_jobs(project_name))}


def cancel_project(project_name: str) -> Dict[str, Any]:
    app = _application()
    manager = _new_manager(app)
    cancelled = manager.cancel_project_jobs(project_name)
    return {"project_name": project_name, **cancelled, "queue": _queue_counts(manager.get_project_jobs(project_name))}


def export_project_results(project_name: str, output_path: os.PathLike[str] | str, output_format: str = "json") -> Dict[str, Any]:
    app = _application()
    data = app.load_project_data(project_name, app.RESULTS_DIR)
    if not data:
        raise JobSpecError(f"No result project found for {project_name!r}.")
    results = data.get("results", [])
    destination = Path(output_path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if output_format == "json":
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(data, handle, indent=2, default=str)
    elif output_format == "csv":
        columns: List[str] = []
        for result in results:
            if isinstance(result, dict):
                for key in result:
                    if key not in columns:
                        columns.append(key)
        with destination.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            for result in results:
                if isinstance(result, dict):
                    writer.writerow({
                        key: json.dumps(value, default=str) if isinstance(value, (dict, list, tuple)) else value
                        for key, value in result.items()
                    })
    else:
        raise JobSpecError("format must be either json or csv.")
    return {"project_name": project_name, "results": len(results), "format": output_format, "output": str(destination)}


def _emit(report: Any) -> None:
    print(json.dumps(report, indent=2, default=str))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Submit and manage BoltzOmics screening jobs from YAML/JSON files.")
    commands = parser.add_subparsers(dest="command", required=True)

    submit = commands.add_parser("submit", help="Validate a manifest and enqueue its jobs without running them.")
    submit.add_argument("manifest")

    run = commands.add_parser("run", help="Submit a manifest, drain the CLI queue, then exit.")
    run.add_argument("manifest")
    run.add_argument("--workers", type=int, help="Number of concurrent BoltzOmics workers (default: manifest workers).")

    worker = commands.add_parser("worker", help="Process the queued CLI jobs, then exit.")
    worker.add_argument("--workers", type=int, help="Concurrent workers (default: detected GPUs recorded by queued jobs).")
    worker.add_argument("--poll-seconds", type=float, default=1.0)

    status = commands.add_parser("status", help="Print persisted queue state as JSON.")
    status.add_argument("--project")

    wait = commands.add_parser("wait", help="Wait until a project's CLI queue has no pending or running jobs.")
    wait.add_argument("--project", required=True)
    wait.add_argument("--timeout", type=float)
    wait.add_argument("--poll-seconds", type=float, default=2.0)

    retry = commands.add_parser("retry", help="Requeue failed jobs for a project.")
    retry.add_argument("--project", required=True)

    cancel = commands.add_parser("cancel", help="Remove pending jobs for a project.")
    cancel.add_argument("--project", required=True)

    export = commands.add_parser("results", help="Export a project's results as JSON or CSV.")
    export.add_argument("--project", required=True)
    export.add_argument("--format", choices=("json", "csv"), default="json")
    export.add_argument("--output", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "submit":
            _emit(submit_spec(args.manifest))
            return 0
        if args.command == "run":
            report = run_spec(args.manifest, workers=args.workers)
            _emit(report)
            return int(report.get("exit_code", 0))
        if args.command == "worker":
            _emit(worker_loop(args.workers, poll_seconds=args.poll_seconds))
            return 0
        if args.command == "status":
            _emit(status_report(args.project))
            return 0
        if args.command == "wait":
            _emit(wait_for_project(args.project, timeout=args.timeout, poll_seconds=args.poll_seconds))
            return 0
        if args.command == "retry":
            _emit(retry_project(args.project))
            return 0
        if args.command == "cancel":
            _emit(cancel_project(args.project))
            return 0
        if args.command == "results":
            _emit(export_project_results(args.project, args.output, args.format))
            return 0
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    except TimeoutError as exc:
        print(str(exc), file=sys.stderr)
        return 124
    except (JobSpecError, OSError, RuntimeError) as exc:
        print(f"BoltzOmics CLI error: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
