#!/usr/bin/env python3
"""Export selected-K native SEISMIC tables to an explicit downstream JSON interface."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def from_frame(frame, sequence, selected_k):
    """Consume native (relationship, K, cluster) columns; retain fractional counts."""
    columns = {}
    for column in frame.columns:
        if not isinstance(column, tuple) or len(column) != 3:
            raise ValueError("Expected native cluster table columns (relationship, K, cluster)")
        relationship, k, cluster = column
        if int(k) != selected_k:
            continue
        relation = str(relationship).lower()
        if relation in ("mutated", "m", "informative", "n"):
            key = ("mutated" if relation in ("mutated", "m") else "informative", int(cluster))
            if key in columns:
                raise ValueError("Duplicate native cluster column")
            columns[key] = column
    for cluster in range(1, selected_k + 1):
        if any((relation, cluster) not in columns for relation in ("mutated", "informative")):
            raise ValueError(f"Missing native mutation/informative counts for K={selected_k}, cluster={cluster}")
    rows = {}
    for index, row in frame.iterrows():
        position = int(index[0] if isinstance(index, tuple) else index)
        if position in rows or not 1 <= position <= len(sequence):
            raise ValueError("Duplicate/out-of-range native local position")
        if isinstance(index, tuple) and len(index) > 1 and str(index[1]).upper().replace("U", "T") != sequence[position - 1]:
            raise ValueError(f"Native table base disagrees with reference at {position}")
        rows[position] = row
    if not rows:
        raise ValueError("Native position table is empty")
    lo, hi = min(rows), max(rows)
    profiles = []
    for cluster in range(1, selected_k + 1):
        rates, info, mutated = [], [], []
        for pos in range(lo, hi + 1):
            n = m = None
            if pos in rows:
                n = float(rows[pos][columns[("informative", cluster)]])
                m = float(rows[pos][columns[("mutated", cluster)]])
                if not math.isfinite(n) or not math.isfinite(m):
                    n = m = None
                elif n < 0 or m < 0 or m > n + 1e-6:
                    raise ValueError("Invalid native mutation/informative counts")
            rates.append(m / n if n is not None and n > 0 else None)
            info.append(n)
            mutated.append(m)
        profiles.append({"cluster": cluster, "selected_k": selected_k, "local_start": lo,
                         "local_end": hi, "sequence": sequence[lo - 1:hi],
                         "contexts": [sequence[max(0, p - 2):min(len(sequence), p + 1)] for p in range(lo, hi + 1)],
                         "mutation_rate": rates,
                         "informative": info, "mutated": mutated})
    return profiles


def export_profiles(summary_file, data_dir, output, loader=None):
    if loader is None:
        from seismicrna.cluster.data import ClusterPositionTableLoader
        loader = ClusterPositionTableLoader
    summary_file, data_dir = Path(summary_file), Path(data_dir).resolve()
    summary = json.loads(summary_file.read_text())
    profiles, statuses = [], []
    for region in summary["regions"]:
        reports = region.get("reports", [])
        if not reports:
            statuses.append({"region": region["name"], "status": region["status"]})
            continue
        reference = data_dir / "references" / (region["name"] + ".fa")
        sequence = "".join(reference.read_text().splitlines()[1:]).upper().replace("U", "T")
        if len(sequence) != region["end"] - region["start"]:
            raise ValueError("Saved search reference length disagrees with genomic mapping")
        for relative in reports:
            report_path = (data_dir / relative).resolve()
            if not report_path.is_relative_to(data_dir):
                raise ValueError("Native report path leaves clustering directory")
            report = json.loads(report_path.read_text())
            k = report["Best number of clusters"]
            if not isinstance(k, int) or k < 1:
                statuses.append({"region": region["name"], "report": relative, "status": "no_valid_model"})
                continue
            table_path = report_path.with_name("cluster-position-table.csv")
            table = loader(table_path)
            for profile in from_frame(table.data, sequence, k):
                # Domain name is native; digest avoids collisions after sanitization.
                identity = hashlib.sha256((relative + str(profile["cluster"])).encode()).hexdigest()[:16]
                lo, hi = profile["local_start"], profile["local_end"]
                if region["strand"] == "+":
                    start, end = region["start"] + lo - 1, region["start"] + hi
                else:
                    start, end = region["end"] - hi, region["end"] - lo + 1
                profile.update(id=f"{region['name']}_k{k}_c{profile['cluster']}_{identity}",
                               sample=summary["sample"], search_region=region["name"], domain=report["Region"],
                               chrom=region["chrom"], start=start, end=end, strand=region["strand"],
                               native_report=relative, native_position_table=str(table_path.relative_to(data_dir)),
                               native_abundance_table=str(report_path.with_name("cluster-abundance-table.csv").relative_to(data_dir)))
                profiles.append(profile)
    payload = {"schema": "insight-cluster-profiles-v1", "sample": summary["sample"],
               "engine": summary["engine"], "source_summary": str(summary_file.resolve()),
               "counts": "native fractional cluster counts; rate=Mutated/Informative",
               "coordinates": "genomic BED0; arrays in transcript 5-prime to 3-prime orientation",
               "profiles": profiles, "statuses": statuses}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("summary", "data-dir", "output"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    export_profiles(args.summary, args.data_dir, args.output)
