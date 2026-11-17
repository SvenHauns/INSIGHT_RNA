"""Configuration and sample-table handling for genomic DMS preprocessing."""
import csv
from pathlib import Path
import re

DEFAULTS = {
    "output": "Preprocessing", "threads": 8, "mem_mb": 64000,
    "min_length": 25, "quality_window": 6, "quality_mean": 25,
    "min_mapq": 25, "max_mismatch_fraction": 0.15,
    "min_score_fraction": 0.3, "min_match_fraction": 0.3,
    "multimap_max": 20, "sjdb_overhang": 149, "sort_mem_mb": 1024,
    "adapter_r1": "", "adapter_r2": "", "poly_g": "auto",
}


def validate_config(config):
    settings = {**DEFAULTS, **config}
    for field in ("samples", "genome_fasta", "annotation_gtf"):
        if not settings.get(field):
            raise ValueError(f"Missing preprocessing configuration: {field}")
        settings[field] = str(Path(settings[field]).expanduser().resolve())
        if not Path(settings[field]).is_file():
            raise ValueError(f"{field} does not exist: {settings[field]}")
    for field in ("output", "star_index", "scratch"):
        if settings.get(field):
            settings[field] = str(Path(settings[field]).expanduser().resolve())
    if settings.get("star_index") and not (Path(settings["star_index"]) / "Genome").is_file():
        raise ValueError("star_index must be a completed STAR index directory containing Genome")
    for field in ("threads", "mem_mb", "min_length", "quality_window", "quality_mean", "multimap_max", "sjdb_overhang", "sort_mem_mb"):
        value = settings[field]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{field} must be a positive integer")
    if not isinstance(settings["min_mapq"], int) or not 0 <= settings["min_mapq"] <= 255:
        raise ValueError("min_mapq must be between 0 and 255")
    for field in ("max_mismatch_fraction", "min_score_fraction", "min_match_fraction"):
        if isinstance(settings[field], bool) or not isinstance(settings[field], (int, float)) or not 0 < settings[field] <= 1:
            raise ValueError(f"{field} must be in (0, 1]")
    if settings["poly_g"] not in ("auto", "yes", "no"):
        raise ValueError("poly_g must be auto, yes or no")
    for field in ("adapter_r1", "adapter_r2"):
        if settings[field] and not re.fullmatch("[ACGTNacgtn]+", settings[field]):
            raise ValueError(f"{field} must be a DNA sequence")
    if settings["threads"] < 2:
        raise ValueError("threads must be at least 2 for streamed alignment and sorting")
    return settings


def load_samples(filename):
    filename = Path(filename)
    result = {}
    paths = set()
    with filename.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"sample", "run", "r1", "r2"}.issubset(reader.fieldnames or []):
            raise ValueError("Samples TSV must contain sample, run, r1, r2 columns")
        for row in reader:
            if not any(row.values()):
                continue
            for field in ("sample", "run"):
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row.get(field, "")):
                    raise ValueError(f"Invalid {field} ID: {row.get(field)!r}")
            if row["run"] in result:
                raise ValueError(f"Duplicate run: {row['run']}")
            for field in ("r1", "r2"):
                if not row.get(field):
                    raise ValueError(f"Missing {field} for {row['run']}")
                path = Path(row[field]).expanduser()
                path = path if path.is_absolute() else filename.resolve().parent / path
                path = path.resolve()
                if re.search(r"[\s,]", str(path)):
                    raise ValueError("STAR FASTQ paths must not contain whitespace or commas")
                if not str(path).endswith((".fastq.gz", ".fq.gz")):
                    raise ValueError(f"Expected gzip FASTQ: {path}")
                if path in paths:
                    raise ValueError(f"FASTQ reused in multiple sample/run entries: {path}")
                paths.add(path)
                row[field] = str(path)
            result[row["run"]] = row
    if not result:
        raise ValueError("Samples TSV contains no runs")
    return result


def samples_to_runs(runs):
    samples = {}
    for run, row in runs.items():
        samples.setdefault(row["sample"], []).append(run)
    return samples
