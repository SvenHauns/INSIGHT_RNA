#!/usr/bin/env python3
"""Create a paired-end samples TSV from an ENA downloader manifest."""
import argparse
import csv
import json
from pathlib import Path
import re

# Verified run-to-biological-replicate mapping from the downloaded GEO metadata.
HEK293 = {
    "SRR34020750": "hek293_5utr_rep2", "SRR34020751": "hek293_5utr_rep2",
    "SRR34020752": "hek293_5utr_rep1", "SRR34020753": "hek293_5utr_rep1",
}


def make_samples(manifest, output, hek293=False):
    manifest, output = Path(manifest).resolve(), Path(output)
    rows = json.loads(manifest.read_text())
    grouped = {}
    for row in rows:
        run = row["run"]
        path = Path(row["path"])
        path = path if path.is_absolute() else manifest.parent / path
        match = re.search(r"_([12])\.fastq\.gz$", path.name)
        if not match:
            raise ValueError(f"Expected paired-end _1/_2.fastq.gz file: {path}")
        mates = grouped.setdefault(run, {})
        if match[1] in mates:
            raise ValueError(f"Duplicate mate {match[1]} for {run}")
        mates[match[1]] = str(path.resolve())
    if not grouped:
        raise ValueError("Manifest contains no runs")
    records = []
    for run, mates in sorted(grouped.items()):
        if set(mates) != {"1", "2"}:
            raise ValueError(f"Run {run} must have both mates")
        if hek293 and run not in HEK293:
            raise ValueError(f"Run {run} is not one of the four standard-condition HEK293 runs")
        records.append({"sample": HEK293[run] if hek293 else run, "run": run, "r1": mates["1"], "r2": mates["2"]})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=["sample", "run", "r1", "r2"], delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hek293-5utr", action="store_true", help="Group the four known technical runs into two biological replicates")
    args = parser.parse_args()
    try:
        records = make_samples(args.manifest, args.output, args.hek293_5utr)
        print(f"Wrote {len(records)} runs across {len({row['sample'] for row in records})} samples to {args.output}")
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"Error: {exc}\n")
