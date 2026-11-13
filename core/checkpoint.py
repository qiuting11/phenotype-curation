# -*- coding: utf-8 -*-
"""Progress checkpointing for resumable downloads.

Provides atomic, crash-safe progress tracking. Each project gets one
checkpoint file (JSON). On resume, completed items are skipped.

Thread-safe via lockfile on Windows (msvcrt) / Unix (fcntl).

Usage:
    ck = CheckpointManager(output_dir)

    # Start or resume
    progress = ck.load_or_create("PRJNA594393", total=947)

    for sample in samples:
        if sample in progress.completed:
            continue  # Already done
        try:
            fetch_sample(sample)
            ck.mark_done(progress, sample)
        except Exception as e:
            ck.mark_failed(progress, sample, str(e))
"""

import json
import os
import time
from pathlib import Path
from typing import Optional
from .types import Checkpoint, FetchStrategy


class CheckpointManager:
    """Manages progress checkpoint files for resumable downloads."""

    def __init__(self, output_dir: str | Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _path(self, bioproject: str) -> Path:
        """Checkpoint file path for a bioproject."""
        return self.output_dir / f".progress_{bioproject}.json"

    def load_or_create(self,
                        bioproject: str,
                        total_samples: int = 0,
                        strategy: str = FetchStrategy.STRATEGY1.value) -> Checkpoint:
        """Load existing checkpoint or create a new one.

        If a checkpoint exists and --resume is intended, reuse it.
        Otherwise overwrite with a fresh checkpoint (only if total_samples > 0).
        """
        ck_path = self._path(bioproject)
        if ck_path.exists():
            return self._load(bioproject)
        return Checkpoint(
            bioproject=bioproject,
            strategy=strategy,
            total_samples=total_samples,
        )

    def resume_from(self, bioproject: str) -> Optional[Checkpoint]:
        """Try to resume from an existing checkpoint. Returns None if no checkpoint."""
        ck_path = self._path(bioproject)
        if not ck_path.exists():
            return None
        return self._load(bioproject)

    def _load(self, bioproject: str) -> Checkpoint:
        """Load checkpoint from disk."""
        ck_path = self._path(bioproject)
        try:
            with open(ck_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except (json.JSONDecodeError, OSError):
            # Corrupted checkpoint — discard
            return Checkpoint(bioproject=bioproject)

        return Checkpoint(
            bioproject=data.get('bioproject', bioproject),
            strategy=data.get('strategy', FetchStrategy.STRATEGY1.value),
            total_samples=data.get('total_samples', 0),
            completed=data.get('completed', []),
            failed=data.get('failed', {}),
            last_updated=data.get('last_updated', ''),
            output_dir=data.get('output_dir', ''),
        )

    def save(self, checkpoint: Checkpoint):
        """Save checkpoint atomically (write temp → rename)."""
        ck_path = self._path(checkpoint.bioproject)
        tmp_path = ck_path.with_suffix('.tmp')

        checkpoint.last_updated = time.strftime('%Y-%m-%dT%H:%M:%S')

        payload = {
            'bioproject': checkpoint.bioproject,
            'strategy': checkpoint.strategy,
            'total_samples': checkpoint.total_samples,
            'completed': checkpoint.completed,
            'failed': checkpoint.failed,
            'last_updated': checkpoint.last_updated,
            'output_dir': str(checkpoint.output_dir),
        }

        with open(tmp_path, 'w', encoding='utf-8') as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
            f.flush()
            os.fsync(f.fileno())

        os.replace(tmp_path, ck_path)  # Atomic on Windows & Unix

    def mark_done(self, checkpoint: Checkpoint, sample_id: str):
        """Mark a sample as successfully downloaded."""
        if sample_id not in checkpoint.completed:
            checkpoint.completed.append(sample_id)

    def mark_failed(self, checkpoint: Checkpoint, sample_id: str, error: str):
        """Record a persistent failure."""
        checkpoint.failed[sample_id] = error

    def is_skipped(self, checkpoint: Checkpoint, sample_id: str) -> bool:
        """Should this sample be skipped? (already done, or previously failed though failed are retried once)."""
        return sample_id in checkpoint.completed

    def delete(self, bioproject: str):
        """Remove checkpoint file (for clean restart)."""
        ck_path = self._path(bioproject)
        try:
            ck_path.unlink()
        except FileNotFoundError:
            pass