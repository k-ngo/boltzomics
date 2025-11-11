import atexit
import json
import os
import threading
import time
import uuid
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Tuple, Callable


@dataclass
class ScreeningJob:
    """Represents a single screening prediction job."""

    job_id: str
    project_name: str
    protein_name: str
    protein_sequence: str
    drug_name: str
    smiles: str
    structure_only: bool
    parameters: Dict[str, Any]
    workspace_name: str
    design_name: str
    status: str = "pending"
    retries: int = 0
    max_attempts: int = 1
    error: Optional[str] = None
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    completed_at: Optional[float] = None
    result: Optional[Dict[str, Any]] = None
    result_metadata: Optional[Dict[str, Any]] = None
    result_committed: bool = False

    def __post_init__(self) -> None:
        self._signature_cache: Optional[str] = None

    @staticmethod
    def _normalize_config(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: ScreeningJob._normalize_config(value[key]) for key in sorted(value)}
        if isinstance(value, (list, tuple)):
            return [ScreeningJob._normalize_config(item) for item in value]
        if isinstance(value, set):
            normalized = [ScreeningJob._normalize_config(item) for item in value]
            return sorted(normalized, key=lambda item: json.dumps(item, sort_keys=True))
        return value

    @classmethod
    def build_signature(
        cls,
        project_name: str,
        protein_name: str,
        protein_sequence: str,
        drug_name: str,
        smiles: str,
        structure_only: bool,
        parameters: Dict[str, Any],
    ) -> str:
        payload = {
            "project_name": project_name,
            "protein_name": protein_name,
            "protein_sequence": protein_sequence,
            "drug_name": drug_name,
            "smiles": smiles,
            "structure_only": structure_only,
            "parameters": cls._normalize_config(parameters or {}),
        }
        normalized = cls._normalize_config(payload)
        return json.dumps(normalized, sort_keys=True, separators=(",", ":"), ensure_ascii=False)

    def configuration_signature(self) -> str:
        if self._signature_cache is None:
            self._signature_cache = self.build_signature(
                project_name=self.project_name,
                protein_name=self.protein_name,
                protein_sequence=self.protein_sequence,
                drug_name=self.drug_name,
                smiles=self.smiles,
                structure_only=self.structure_only,
                parameters=self.parameters,
            )
        return self._signature_cache

    def set_cached_signature(self, signature: str) -> None:
        self._signature_cache = signature

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ScreeningJob":
        return cls(**data)


class ScreeningJobManager:
    """Simple background job manager for screening predictions."""

    def __init__(self, state_dir: str):
        self.state_dir = state_dir
        os.makedirs(self.state_dir, exist_ok=True)
        self.state_path = os.path.join(self.state_dir, "job_state.json")
        self.jobs: Dict[str, ScreeningJob] = {}
        self.lock = threading.Lock()
        self.new_job_event = threading.Event()
        self.shutdown_event = threading.Event()
        self.worker_thread: Optional[threading.Thread] = None
        self.executor: Optional[Callable[[ScreeningJob], Tuple[Dict[str, Any], Dict[str, Any]]]] = None
        self._load_state()
        self._atexit_registered = False
        self._register_atexit()
        self._last_enqueue_report: Dict[str, List[str]] = {"enqueued": [], "skipped": []}

    def _register_atexit(self) -> None:
        if self._atexit_registered:
            return

        def _cleanup() -> None:
            try:
                self.shutdown(purge_pending=True)
            except Exception:
                pass

        atexit.register(_cleanup)
        self._atexit_registered = True

    def register_executor(
        self,
        executor: Callable[[ScreeningJob], Tuple[Dict[str, Any], Dict[str, Any]]],
    ) -> None:
        """Register the callable used to execute jobs."""
        self.executor = executor
        self._ensure_worker()

    def _find_jobs_by_signature_locked(self, signature: str) -> Tuple[Optional[ScreeningJob], List[str]]:
        active_job: Optional[ScreeningJob] = None
        redundant_ids: List[str] = []
        for job_id, existing in self.jobs.items():
            if existing.configuration_signature() != signature:
                continue
            if existing.status in {"failed", "cancelled"}:
                redundant_ids.append(job_id)
            else:
                active_job = existing
        return active_job, redundant_ids

    def enqueue_jobs(self, jobs: Iterable[ScreeningJob]) -> List[ScreeningJob]:
        """Add jobs to the queue."""
        enqueued: List[ScreeningJob] = []
        skipped_signatures: List[str] = []
        with self.lock:
            updated = False
            for job in jobs:
                signature = job.configuration_signature()
                active_job, redundant_ids = self._find_jobs_by_signature_locked(signature)
                if active_job:
                    skipped_signatures.append(signature)
                    continue
                for redundant_id in redundant_ids:
                    if redundant_id in self.jobs:
                        del self.jobs[redundant_id]
                        updated = True
                if job.job_id in self.jobs:
                    del self.jobs[job.job_id]
                    updated = True
                self.jobs[job.job_id] = job
                enqueued.append(job)
                updated = True
            if updated:
                self._persist_state_locked()
            # Store a lightweight summary for diagnostic UI
            deduped_skipped = list(dict.fromkeys(skipped_signatures))
            self._last_enqueue_report = {
                "enqueued": [job.job_id for job in enqueued],
                "skipped": deduped_skipped,
            }
        if enqueued:
            self.new_job_event.set()
            self._ensure_worker()
        return enqueued

    def get_project_jobs(self, project_name: str) -> List[ScreeningJob]:
        with self.lock:
            return [job for job in self.jobs.values() if job.project_name == project_name]

    def get_project_summary(self, project_name: str) -> Dict[str, Any]:
        jobs = self.get_project_jobs(project_name)
        summary = {
            "pending": len([j for j in jobs if j.status == "pending"]),
            "running": len([j for j in jobs if j.status == "running"]),
            "success": len([j for j in jobs if j.status == "success"]),
            "failed": len([j for j in jobs if j.status == "failed"]),
            "cancelled": len([j for j in jobs if j.status == "cancelled"]),
            "total": len(jobs),
            "active_job": next((j for j in jobs if j.status == "running"), None),
        }

        if not jobs:
            summary.update(
                {
                    "elapsed_seconds": None,
                    "eta_seconds": None,
                    "average_duration_seconds": None,
                }
            )
            return summary

        now = time.time()

        # Compute total elapsed time by summing up individual job elapsed times
        elapsed_seconds: Optional[float] = None
        total_elapsed = 0.0
        has_elapsed_data = False

        for job in jobs:
            if job.completed_at and job.started_at and job.completed_at >= job.started_at:
                # Completed jobs: use actual duration
                total_elapsed += (job.completed_at - job.started_at)
                has_elapsed_data = True
            elif job.status == "running" and job.started_at:
                # Running jobs: use time elapsed so far
                total_elapsed += max(0.0, now - job.started_at)
                has_elapsed_data = True
            elif job.status == "failed" and job.started_at:
                # Failed jobs: count time from start to completion (if available) or to now
                if job.completed_at:
                    total_elapsed += max(0.0, job.completed_at - job.started_at)
                else:
                    total_elapsed += max(0.0, now - job.started_at)
                has_elapsed_data = True

        if has_elapsed_data:
            elapsed_seconds = total_elapsed

        durations = [
            job.completed_at - job.started_at
            for job in jobs
            if job.completed_at and job.started_at and job.completed_at >= job.started_at
        ]
        average_duration = sum(durations) / len(durations) if durations else None

        eta_seconds: Optional[float] = None
        if average_duration and average_duration > 0:
            eta_seconds = 0.0
            running_jobs = [job for job in jobs if job.status == "running" and job.started_at]
            for job in running_jobs:
                elapsed_running = max(0.0, now - job.started_at)
                remaining = max(0.0, average_duration - elapsed_running)
                eta_seconds += remaining
            eta_seconds += average_duration * summary["pending"]
            if eta_seconds == 0.0:
                eta_seconds = None

        summary.update(
            {
                "elapsed_seconds": elapsed_seconds,
                "eta_seconds": eta_seconds,
                "average_duration_seconds": average_duration,
            }
        )
        return summary

    def get_uncommitted_results(self, project_name: str) -> List[Tuple[ScreeningJob, Dict[str, Any], Dict[str, Any]]]:
        jobs = self.get_project_jobs(project_name)
        ready: List[Tuple[ScreeningJob, Dict[str, Any], Dict[str, Any]]] = []
        for job in jobs:
            if job.status == "success" and job.result and not job.result_committed:
                ready.append((job, job.result, job.result_metadata or {}))
        return ready

    def mark_results_committed(self, job_ids: Iterable[str]) -> None:
        with self.lock:
            updated = False
            for job_id in job_ids:
                job = self.jobs.get(job_id)
                if job and not job.result_committed:
                    job.result_committed = True
                    updated = True
            if updated:
                self._persist_state_locked()

    def has_job_with_signature(self, signature: str, include_failed: bool = False) -> Optional[ScreeningJob]:
        with self.lock:
            for job in self.jobs.values():
                if job.configuration_signature() != signature:
                    continue
                if include_failed or job.status not in {"failed", "cancelled"}:
                    return job
        return None

    def _purge_pending_jobs_locked(self, project_name: Optional[str] = None) -> List[str]:
        removed_ids: List[str] = []
        for job_id, job in list(self.jobs.items()):
            if job.status != "pending":
                continue
            if project_name and job.project_name != project_name:
                continue
            removed_ids.append(job_id)
        for job_id in removed_ids:
            del self.jobs[job_id]
        return removed_ids

    def cancel_project_jobs(self, project_name: str) -> Dict[str, int]:
        with self.lock:
            removed_ids = self._purge_pending_jobs_locked(project_name)
            if removed_ids:
                self._persist_state_locked()
        if removed_ids:
            self.new_job_event.set()
        return {
            "pending_removed": len(removed_ids),
            "running_cancelled": 0,
        }

    def shutdown(self, purge_pending: bool = False) -> None:
        if purge_pending:
            with self.lock:
                removed_ids = self._purge_pending_jobs_locked()
                if removed_ids:
                    self._persist_state_locked()
        self.shutdown_event.set()
        self.new_job_event.set()
        if self.worker_thread and self.worker_thread.is_alive():
            self.worker_thread.join(timeout=5)

    def _ensure_worker(self) -> None:
        if self.worker_thread and self.worker_thread.is_alive():
            return
        if self.executor is None:
            return
        self.worker_thread = threading.Thread(target=self._worker_loop, daemon=True)
        self.worker_thread.start()

    def _worker_loop(self) -> None:
        while not self.shutdown_event.is_set():
            job = self._start_next_job()
            if job is None:
                self.new_job_event.wait(timeout=1.0)
                self.new_job_event.clear()
                continue
            try:
                result, metadata = self.executor(job)
            except Exception as exc:
                with self.lock:
                    job.retries += 1
                    job.status = "failed"
                    job.error = str(exc)
                    job.completed_at = time.time()
                    self._persist_state_locked()
                self.new_job_event.set()
                continue
            with self.lock:
                job.status = "success"
                job.error = None
                job.completed_at = time.time()
                job.result = result
                job.result_metadata = metadata
                self._persist_state_locked()
            self.new_job_event.set()

    def _start_next_job(self) -> Optional[ScreeningJob]:
        with self.lock:
            for job in self.jobs.values():
                if job.status == "pending":
                    job.status = "running"
                    job.started_at = time.time()
                    self._persist_state_locked()
                    return job
        return None

    def _load_state(self) -> None:
        if not os.path.exists(self.state_path):
            return
        try:
            with open(self.state_path, "r") as f:
                data = json.load(f)
        except Exception:
            return
        jobs = data.get("jobs", [])
        for entry in jobs:
            try:
                job = ScreeningJob.from_dict(entry)
                self.jobs[job.job_id] = job
            except Exception:
                continue

    def _persist_state_locked(self) -> None:
        state = {
            "jobs": [job.to_dict() for job in self.jobs.values()],
            "updated_at": time.time(),
            "version": 1,
        }
        tmp_path = os.path.join(self.state_dir, f"job_state_{uuid.uuid4().hex}.tmp")
        with open(tmp_path, "w") as f:
            json.dump(state, f, indent=2)
        os.replace(tmp_path, self.state_path)


__all__ = ["ScreeningJob", "ScreeningJobManager"]
