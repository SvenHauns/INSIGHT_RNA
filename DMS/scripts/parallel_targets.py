"""Bounded, ordered target execution with a process-local BAM handle."""
from concurrent.futures import ProcessPoolExecutor
from collections import deque
from contextlib import contextmanager
import atexit
import multiprocessing
import os
from pathlib import Path
import tempfile

_BAM = None


def close_bam():
    global _BAM
    if _BAM is not None:
        _BAM.close()
        _BAM = None


def initialize_bam(path):
    global _BAM
    import pysam
    close_bam()
    try:
        _BAM = pysam.AlignmentFile(path)
        _BAM.check_index()
    except Exception as exc:
        close_bam()
        raise RuntimeError(f"Cannot open indexed BAM {path}: {exc}") from exc
    atexit.register(close_bam)


def get_bam():
    if _BAM is None:
        raise RuntimeError("Target worker BAM was not initialized")
    return _BAM


def ordered_results(function, tasks, workers, initializer=None, initargs=()):
    """Keep at most 2*workers results in flight, including slow early targets."""
    if workers < 1:
        raise ValueError("workers must be positive")
    if workers == 1:
        if initializer:
            initializer(*initargs)
        try:
            for task in tasks:
                yield function(task)
        finally:
            close_bam()
        return
    with ProcessPoolExecutor(max_workers=workers,
                             mp_context=multiprocessing.get_context("spawn"),
                             initializer=initializer, initargs=initargs) as pool:
        pending = deque()
        tasks = iter(tasks)
        try:
            for _ in range(2 * workers):
                task = next(tasks, None)
                if task is None:
                    break
                pending.append(pool.submit(function, task))
            while pending:
                yield pending.popleft().result()
                task = next(tasks, None)
                if task is not None:
                    pending.append(pool.submit(function, task))
        finally:
            for future in pending:
                future.cancel()


@contextmanager
def atomic_output(path, mode="w"):
    """Never append to an old or interrupted aggregate output."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(descriptor, mode) as stream:
            yield stream
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
