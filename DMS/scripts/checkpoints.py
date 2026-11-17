"""Atomic target checkpoints outside Snakemake's replaceable rule outputs."""
import hashlib
import json
from pathlib import Path
import pickle

from parallel_targets import atomic_output, ordered_results


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    identity = {'path': str(path), 'size': stat.st_size}
    # Regenerated BEDs/coverage lists/FASTA slices can have identical content but
    # new timestamps. Content hashes let their dependent profile caches survive.
    if path.suffix in {'.bed', '.txt', '.json', '.fai'}:
        hasher = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                hasher.update(block)
        identity['sha256'] = hasher.hexdigest()
    else:
        # Never reread a large genomic FASTA or BAM merely to resume a job.
        identity['mtime_ns'] = stat.st_mtime_ns
    return identity


def checkpoint_signature(stage, inputs, settings=None):
    scripts = Path(__file__).resolve().parent
    code = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(scripts.glob('*.py'))}
    return digest({'schema': 1, 'stage': stage, 'inputs': [file_identity(path) for path in inputs],
                   'settings': settings, 'code': code})


def read_result(path, key):
    try:
        with path.open('rb') as stream:
            payload = pickle.load(stream)
        if payload['key'] == key:
            return True, payload['result']
    except (OSError, EOFError, pickle.UnpicklingError, KeyError, TypeError, ValueError):
        pass
    return False, None


def checkpointed_results(function, tasks, workers, initializer, initargs,
                         stage, inputs, checkpoint_dir='.dms-checkpoints', resume=True):
    if not resume:
        yield from ordered_results(function, tasks, workers, initializer, initargs)
        return
    directory = Path(checkpoint_dir) / stage / checkpoint_signature(stage, inputs)
    plan = []
    for task in tasks:
        key = digest(task)
        path = directory / key[:2] / (key + '.pkl')
        valid, _ = read_result(path, key)
        plan.append((task, key, path, valid))
    reused = sum(valid for _, _, _, valid in plan)
    print(f'{stage}: reusing {reused}/{len(plan)} completed target checkpoints', flush=True)
    missing = (task for task, _, _, valid in plan if not valid)
    results = iter(ordered_results(function, missing, workers, initializer, initargs))
    try:
        for task, key, path, valid in plan:
            if valid:
                # The same workflow's Snakemake lock prevents concurrent writers.
                valid, result = read_result(path, key)
                if not valid:
                    raise RuntimeError(f'Checkpoint changed while reading: {path}; resubmit after checking concurrent jobs')
            else:
                result = next(results)
                with atomic_output(path, 'wb') as stream:
                    pickle.dump({'key': key, 'result': result}, stream)
            yield result
    finally:
        results.close()
