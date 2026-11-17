#!/usr/bin/env python3
"""Optional native SEISMIC-RNA scanning/clustering of strand-aware genomic regions."""
import argparse
from dataclasses import dataclass
import json
import os
import shutil
import uuid
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

from checkpoints import checkpoint_signature
from parallel_targets import atomic_output


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str
    strand: str


def enabled(value):
    if isinstance(value, bool):
        return value
    text = str(value).lower()
    if text not in {"true", "false", "1", "0", "yes", "no"}:
        raise ValueError("cluster_dms must be true or false")
    return text in {"true", "1", "yes"}


def read_regions(path, max_length=None):
    regions = []
    names = set()
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 6:
                raise ValueError(f"BED line {number}: expected BED6 (chrom, start, end, name, score, strand)")
            chrom, start, end, name, _, strand = fields
            start, end = int(start), int(end)
            if start < 0 or end <= start or (max_length is not None and end - start > max_length):
                raise ValueError(f"BED line {number}: search region must have positive length and nonnegative coordinates (limit: {max_length})")
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or name in names:
                raise ValueError(f"Invalid or duplicate region ID: {name}")
            if strand not in {"+", "-"} or re.search(r"[\s:]", chrom):
                raise ValueError(f"BED line {number}: explicit strand (+/-) and a valid chromosome name are required")
            names.add(name)
            regions.append(Region(chrom, start, end, name, strand))
    if not regions:
        raise ValueError("No clustering regions found")
    return regions


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def crop_alignment(fields, region):
    """Clip a SAM alignment to a contiguous window, preserving SEQ and qualities.

    Coordinates become window-local. Reverse-strand windows also reverse CIGAR,
    SEQ, QUAL and orientation so DMS A/C selection is on the transcript strand.
    MD/NM and mate coordinates must be regenerated after this operation.
    """
    flag = int(fields[1])
    if flag & 0xB04 or fields[2] != region.chrom or fields[5] == "*" or fields[9] == "*":
        return None
    operations = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", fields[5])]
    if "".join(f"{n}{op}" for n, op in operations) != fields[5]:
        raise ValueError(f"Unsupported CIGAR: {fields[5]}")
    ref, query = int(fields[3]) - 1, 0
    selected = []
    first_query = last_query = first_ref = last_ref = None
    for length, op in operations:
        if op in "M=X":
            lo, hi = max(ref, region.start), min(ref + length, region.end)
            if lo < hi:
                qlo, qhi = query + lo - ref, query + hi - ref
                if first_query is None:
                    first_query, first_ref = qlo, lo
                last_query, last_ref = qhi, hi
                selected.append((hi - lo, op))
        elif op in "DN":
            lo, hi = max(ref, region.start), min(ref + length, region.end)
            if lo < hi:
                selected.append((hi - lo, op))
        elif op == "I" and region.start <= ref < region.end:
            selected.append((length, op))
        if op in "MDN=X":
            ref += length
        if op in "MIS=X":
            query += length
    if first_query is None:
        return None
    # Non-matching operations at boundaries cannot anchor a local alignment.
    while selected and selected[0][1] not in "M=X":
        selected.pop(0)
    while selected and selected[-1][1] not in "M=X":
        selected.pop()
    cigar = []
    if operations[0][1] == "H":
        cigar.append(operations[0])
    if first_query:
        cigar.append((first_query, "S"))
    cigar += selected
    if last_query < query:
        cigar.append((query - last_query, "S"))
    if operations[-1][1] == "H":
        cigar.append(operations[-1])
    result = list(fields)
    result[2] = region.name
    result[3] = str(first_ref - region.start + 1)
    result[6:9] = ["*", "0", "0"]
    # fixmate will reconstruct local mate flags/coordinates after name sorting.
    flag &= ~(2 | 8 | 32)
    if region.strand == "-":
        cigar.reverse()
        result[3] = str(region.end - last_ref + 1)
        result[9] = reverse_complement(result[9])
        if result[10] != "*":
            result[10] = result[10][::-1]
        flag ^= 16
    result[1] = str(flag)
    result[5] = "".join(f"{n}{op}" for n, op in cigar)
    # Alignment tags derived from the old coordinates/sequence are invalid now.
    result = result[:11] + [tag for tag in result[11:] if tag.split(":", 1)[0] not in {"MD", "NM", "MC", "SA", "AS", "XS", "RG"}]
    return result


def run(command, log, **kwargs):
    log.write(shlex.join(map(str, command)) + "\n")
    log.flush()
    return subprocess.run(list(map(str, command)), stderr=log, check=True, **kwargs)


def prepare_region(bam, fasta, region, directory, log):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    local_fasta = directory / "reference.fa"
    sequence = subprocess.check_output(["samtools", "faidx", str(fasta),
                                        f"{region.chrom}:{region.start + 1}-{region.end}"], text=True)
    sequence = "".join(sequence.splitlines()[1:]).upper()
    if len(sequence) != region.end - region.start:
        raise ValueError(f"Reference does not cover region {region.name}")
    if region.strand == "-":
        sequence = reverse_complement(sequence)
    local_fasta.write_text(f">{region.name}\n{sequence}\n")
    run(["samtools", "faidx", local_fasta], log, stdout=log)
    local_sam = directory / "cropped.sam"
    count = 0
    with local_sam.open("w") as output:
        output.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{region.name}\tLN:{len(sequence)}\n")
        process = subprocess.Popen(["samtools", "view", str(bam), f"{region.chrom}:{region.start + 1}-{region.end}"],
                                   stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            for line in process.stdout:
                record = crop_alignment(line.rstrip("\n").split("\t"), region)
                if record:
                    output.write("\t".join(record) + "\n")
                    count += 1
            if process.wait():
                raise subprocess.CalledProcessError(process.returncode, process.args)
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            process.wait()
    if not count:
        return local_fasta, None, count
    named, fixed, sorted_bam = [directory / name for name in ("named.bam", "fixed.bam", "sorted.bam")]
    run(["samtools", "sort", "-n", "-o", named, local_sam], log, stdout=log)
    run(["samtools", "fixmate", named, fixed], log, stdout=log)
    run(["samtools", "sort", "-o", sorted_bam, fixed], log, stdout=log)
    local_bam = directory / (region.name + ".bam")
    with local_bam.open("wb") as output:
        run(["samtools", "calmd", "-b", sorted_bam, local_fasta], log, stdout=output)
    run(["samtools", "index", local_bam], log, stdout=log)
    for file in (local_sam, named, fixed, sorted_bam):
        file.unlink()
    return local_fasta, local_bam, count


def seismic_commands(fasta, bam, output, sample, region, cpus, max_clusters, min_info,
                     scan=True, pair_fdr=None, detect_fdr=None, merge_fdr=None, validate_gaps=False):
    # BAM's parent layout determines sample name in SEISMIC-RNA.
    idmut = Path(output) / sample / "idmut" / region
    filtered = Path(output) / sample / "filter" / region
    commands = [
        ["seismic", "idmut", str(fasta), str(bam), "--out-dir", str(output),
         "--min-reads", "1", "--write-read-names", "--num-cpus", str(cpus)],
        ["seismic", "filter", str(idmut), "--probe", "DMS", "--min-ninfo-pos", str(min_info),
         "--num-cpus", str(cpus)],
        ["seismic", "cluster", str(filtered), "--max-clusters", str(max_clusters),
         "--cluster-pos-table", "--cluster-abundance-table", "--num-cpus", str(cpus)],
    ]
    if not scan:
        return commands
    workflow = ["seismic", "wf", str(fasta), str(idmut / "idmut-report.json"), "--out-dir", str(output),
                "--scan", "--cluster", "--no-fold", "--no-draw", "--no-collate",
                "--probe", "DMS", "--min-ninfo-pos", str(min_info),
                "--max-clusters", str(max_clusters), "--cluster-pos-table",
                "--cluster-abundance-table", "--num-cpus", str(cpus)]
    for name, value in (("pair-fdr", pair_fdr), ("detect-fdr", detect_fdr), ("merge-fdr", merge_fdr)):
        if value is not None:
            if not 0 < value <= 1:
                raise ValueError(f"{name} must be in (0, 1]")
            workflow += ["--" + name, str(value)]
    if validate_gaps:
        workflow += ["--validate-gaps"]
    # Let the native workflow manage filterscan -> clusterscan data dependencies,
    # rather than guessing the scan report schema or constructing its domains here.
    return [commands[0], workflow]


def cluster_region(bam, staged, region, sample, output, temporary, log, cpus,
                   max_clusters, min_info, scan, pair_fdr, detect_fdr, merge_fdr,
                   validate_gaps):
    native = output / 'native'
    work = Path(temporary) / region.name / sample / "align"
    local_fasta, local_bam, count = prepare_region(bam, staged, region, work, log)
    entry = {**region.__dict__, "alignments": count, "status": "no_reads", "best_k": None}
    if local_bam is not None:
        # Persist target sequences and genomic coordinate mapping; BAMs remain internal.
        saved_ref = output / "references" / (region.name + ".fa")
        saved_ref.parent.mkdir(exist_ok=True)
        saved_ref.write_text(local_fasta.read_text())
        commands = seismic_commands(saved_ref, local_bam, native, sample, region.name, cpus, max_clusters, min_info,
                                    scan, pair_fdr, detect_fdr, merge_fdr, validate_gaps)
        run(commands[0], log, stdout=log)
        idmut_report = native / sample / "idmut" / region.name / "idmut-report.json"
        idmut_data = json.loads(idmut_report.read_text())
        if idmut_data.get("Number of reads processed successfully") == 0:
            entry["status"] = "no_usable_reads"
            return entry
        run(commands[1], log, stdout=log)
        if scan:
            scan_reports = list((native / sample / "filterscan" / region.name).rglob("*.json"))
            if not scan_reports:
                raise RuntimeError(f"No native filterscan report for {region.name}; inspect commands.log")
            reports = list((native / sample / "cluster" / region.name).rglob("cluster-report.json"))
            entry["scan_reports"] = [str(p.relative_to(output)) for p in scan_reports]
            entry["clusterscan_reports"] = [str(p.relative_to(output))
                                           for p in (native / sample / "clusterscan" / region.name).rglob("*.json")]
            entry["reports"] = [str(p.relative_to(output)) for p in reports]
            entry["best_k"] = [json.loads(p.read_text())["Best number of clusters"] for p in reports]
            entry["status"] = "completed" if reports else "scan_completed_without_cluster_reports"
            if reports and not all(isinstance(k, int) and 1 <= k <= max_clusters for k in entry["best_k"]):
                entry["status"] = "no_valid_model"
            return entry
        filter_reports = list((native / sample / "filter" / region.name).rglob("filter-report.json"))
        if not filter_reports:
            raise RuntimeError(f"No filter report for {region.name}; inspect commands.log")
        filter_data = [json.loads(p.read_text()) for p in filter_reports]
        if all(data.get("Number of positions kept after filtering") == 0 or
               data.get("Number of reads kept after filtering") == 0 for data in filter_data):
            entry["status"] = "insufficient_after_filtering"
            entry["filter_reports"] = [str(p.relative_to(output)) for p in filter_reports]
            return entry
        run(commands[2], log, stdout=log)
        reports = list((native / sample / "cluster" / region.name).rglob("cluster-report.json"))
        if not reports:
            raise RuntimeError(f"No cluster report for {region.name}; inspect {output / 'commands.log'} for filtered/insufficient data or errors")
        entry["status"] = "completed"
        entry["reports"] = [str(p.relative_to(output)) for p in reports]
        entry["best_k"] = [json.loads(p.read_text())["Best number of clusters"] for p in reports]
        if not all(isinstance(k, int) and 1 <= k <= max_clusters for k in entry["best_k"]):
            entry["status"] = "no_valid_model"
    return entry

def region_artifacts(output, sample, name):
    paths = []
    native_sample = output / 'native' / sample
    if native_sample.exists():
        for stage in native_sample.iterdir():
            region_dir = stage / name
            if region_dir.is_dir():
                paths.extend(path for path in region_dir.rglob('*') if path.is_file())
    reference = output / 'references' / (name + '.fa')
    if reference.is_file():
        paths.append(reference)
    return {str(path.relative_to(output)): [path.stat().st_size, path.stat().st_mtime_ns]
            for path in paths}


def completed_region(marker, output, region):
    try:
        payload = json.loads(marker.read_text())
        entry = payload['entry']
        if any(entry[key] != value for key, value in region.__dict__.items()):
            return None
        for relative, expected in payload['artifacts'].items():
            path = (output / relative).resolve()
            if not path.is_relative_to(output.resolve()):
                return None
            stat = path.stat()
            if [stat.st_size, stat.st_mtime_ns] != expected:
                return None
        return entry
    except (OSError, ValueError, KeyError, TypeError):
        return None


def clear_unfinished_region(output, sample, region):
    """Discard only this cache's uncommitted region, never the source BAM."""
    native_sample = output / 'native' / sample
    if native_sample.exists():
        for stage in native_sample.iterdir():
            region_dir = stage / region.name
            if region_dir.is_dir():
                shutil.rmtree(region_dir)
    for suffix in ('.fa', '.fa.fai'):
        (output / 'references' / (region.name + suffix)).unlink(missing_ok=True)


def link_or_copy(source, destination):
    try:
        return os.link(source, destination)
    except OSError:
        return shutil.copy2(source, destination)


def publish_cluster_output(cached, output):
    """Publish hardlinks when possible; checkpoints survive output cleanup."""
    output.parent.mkdir(parents=True, exist_ok=True)
    staged = output.parent / ('.' + output.name + '.publish-' + uuid.uuid4().hex)
    try:
        shutil.copytree(cached, staged, copy_function=link_or_copy)
        if output.exists():
            # Preserve old direct-invocation output instead of mixing native runs.
            archived = output.parent / ('.' + output.name + '.previous-' + uuid.uuid4().hex)
            output.rename(archived)
        staged.rename(output)
    finally:
        if staged.exists():
            shutil.rmtree(staged)


def cluster_sample(bam, fasta, bed, sample, output, summary, cpus=4, max_clusters=3, min_info=1000,
                   scan=True, pair_fdr=None, detect_fdr=None, merge_fdr=None, validate_gaps=False,
                   checkpoint_dir='.dms-checkpoints', resume=True):
    if max_clusters < 1 or cpus < 1 or min_info < 1:
        raise ValueError('cpus, max_clusters and min_info must be positive')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample):
        raise ValueError('Invalid sample ID')
    regions = read_regions(bed)
    output = Path(output).resolve()
    version = subprocess.check_output(['seismic', '--version'], text=True).strip()
    if not re.search(r'\b0\.26\.\d+\b', version):
        raise ValueError(f'Expected SEISMIC-RNA 0.26.x; got {version}')
    parameters = {'sample': sample, 'engine': version, 'max_clusters': max_clusters,
                  'min_info': min_info, 'scan': scan, 'pair_fdr': pair_fdr,
                  'detect_fdr': detect_fdr, 'merge_fdr': merge_fdr, 'validate_gaps': validate_gaps}
    signature = checkpoint_signature('clustering', [bam, str(bam) + '.bai', fasta, bed], parameters)
    cache = Path(checkpoint_dir).resolve() / 'clustering' / sample / signature
    if cache.is_relative_to(output):
        raise ValueError('Clustering checkpoints must be outside the final output directory')
    cached_output = cache / 'data'
    (cached_output / 'native').mkdir(parents=True, exist_ok=True)
    markers = cache / 'regions'
    markers.mkdir(exist_ok=True)
    statuses = []
    # Staged reference/BAM files are disposable; native outputs/checkpoints persist.
    with tempfile.TemporaryDirectory(prefix='insight-cluster-') as temporary:
        staged = Path(temporary) / 'genome.fa'
        staged.symlink_to(Path(fasta).resolve())
        subprocess.run(['samtools', 'faidx', str(staged)], check=True)
        bam_lengths = {}
        header = subprocess.check_output(['samtools', 'view', '-H', str(bam)], text=True)
        for line in header.splitlines():
            if line.startswith('@SQ\t'):
                tags = dict(field.split(':', 1) for field in line.split('\t')[1:])
                bam_lengths[tags['SN']] = int(tags['LN'])
        with Path(str(staged) + '.fai').open() as stream:
            fasta_lengths = {parts[0]: int(parts[1]) for line in stream if (parts := line.split('\t'))}
        for region in regions:
            if region.chrom not in bam_lengths or region.chrom not in fasta_lengths:
                raise ValueError(f'Region {region.name}: chromosome missing from BAM or FASTA')
            if bam_lengths[region.chrom] != fasta_lengths[region.chrom] or region.end > bam_lengths[region.chrom]:
                raise ValueError(f'Region {region.name}: BAM/FASTA length mismatch or coordinates outside reference')
        with (cached_output / 'commands.log').open('a') as log:
            for number, region in enumerate(regions, 1):
                marker = markers / (region.name + '.json')
                entry = completed_region(marker, cached_output, region) if resume else None
                if entry is not None:
                    print(f'{sample}: reusing region {number}/{len(regions)} {region.name}', flush=True)
                else:
                    marker.unlink(missing_ok=True)
                    clear_unfinished_region(cached_output, sample, region)
                    entry = cluster_region(bam, staged, region, sample, cached_output, temporary, log,
                                           cpus, max_clusters, min_info, scan, pair_fdr,
                                           detect_fdr, merge_fdr, validate_gaps)
                    with atomic_output(marker) as stream:
                        json.dump({'entry': entry, 'artifacts': region_artifacts(cached_output, sample, region.name)}, stream)
                    print(f'{sample}: checkpointed region {number}/{len(regions)} {region.name}: {entry["status"]}', flush=True)
                statuses.append(entry)
    publish_cluster_output(cached_output, output)
    with atomic_output(summary) as stream:
        json.dump({**parameters, 'source_bam': str(Path(bam).resolve()),
                   'source_fasta': str(Path(fasta).resolve()), 'regions': statuses}, stream, indent=2)
        stream.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("bam", "fasta", "regions", "sample", "output", "summary"):
        parser.add_argument("--" + field, required=True)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--max-clusters", type=int, default=3)
    parser.add_argument("--min-info", type=int, default=1000)
    parser.add_argument("--scan", action=argparse.BooleanOptionalAction, default=True)
    for name in ("pair-fdr", "detect-fdr", "merge-fdr"):
        parser.add_argument("--" + name, type=float)
    parser.add_argument("--validate-gaps", action="store_true")
    parser.add_argument("--checkpoint-dir", default=".dms-checkpoints")
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    cluster_sample(args.bam, args.fasta, args.regions, args.sample, args.output, args.summary,
                   args.cpus, args.max_clusters, args.min_info, args.scan,
                   args.pair_fdr, args.detect_fdr, args.merge_fdr, args.validate_gaps,
                   args.checkpoint_dir, not args.no_resume)
