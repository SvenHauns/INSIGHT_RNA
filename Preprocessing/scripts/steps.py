"""Small, testable subprocess steps; invoked by Snakemake in its Conda environment."""
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile


def parents(*paths):
    for path in paths:
        Path(path).parent.mkdir(parents=True, exist_ok=True)


def run(command, log, stdout=None):
    log.write((shlex.join(map(str, command)) + "\n").encode())
    log.flush()
    subprocess.run(list(map(str, command)), stdout=stdout or log, stderr=log, check=True)


def pipeline(commands, log, stdout=None):
    """Stream BAMs and fail if any producer or consumer fails, not just the last."""
    processes = []
    previous = None
    try:
        for index, command in enumerate(commands):
            command = list(map(str, command))
            log.write((shlex.join(command) + "\n").encode())
            log.flush()
            process = subprocess.Popen(command, stdin=previous,
                                       stdout=subprocess.PIPE if index < len(commands) - 1 else (stdout or log),
                                       stderr=log)
            if previous is not None:
                previous.close()
            previous = process.stdout
            processes.append(process)
        for process in reversed(processes):
            process.wait()
        for command, process in zip(commands, processes):
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode, command)
    finally:
        if previous is not None:
            previous.close()
        for process in processes:
            if process.poll() is None:
                process.terminate()
        for process in processes:
            process.wait()


def prepare_reference(source, gtf, target, report):
    source, target = Path(source).resolve(), Path(target).absolute()
    if source == target:
        raise ValueError("Use a separate source FASTA; the workflow owns its staged reference path")
    parents(target, report)
    if target.exists() and not target.is_symlink():
        raise ValueError(f"Refusing to replace an existing real FASTA: {target}")
    target.unlink(missing_ok=True)
    target.symlink_to(source)
    subprocess.run(["samtools", "faidx", str(target)], check=True)
    contigs = {}
    with Path(str(target) + ".fai").open() as stream:
        for line in stream:
            fields = line.split("\t")
            contigs[fields[0]] = int(fields[1])
    annotation_contigs = set()
    with Path(gtf).open() as stream:
        for line in stream:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError("Annotation must be an uncompressed 9-column GTF")
                annotation_contigs.add(fields[0])
    if not annotation_contigs:
        raise ValueError("No annotation records found in GTF")
    missing = annotation_contigs - contigs.keys()
    if missing:
        raise ValueError("GTF contigs missing from FASTA (check chr/NC_ names and assembly): " + ", ".join(sorted(missing)[:10]))
    Path(report).write_text(json.dumps({"source_fasta": str(source), "annotation": str(gtf),
                                       "contigs": len(contigs), "annotation_contigs": len(annotation_contigs)}, indent=2) + "\n")


def index_reference(fasta, gtf, output, threads, settings, logfile):
    parents(logfile)
    Path(output).mkdir(parents=True, exist_ok=True)
    with open(logfile, "wb") as log:
        run(["STAR", "--runMode", "genomeGenerate", "--runThreadN", threads,
             "--genomeDir", output, "--genomeFastaFiles", fasta,
             "--sjdbGTFfile", gtf, "--sjdbOverhang", settings["sjdb_overhang"],
             "--limitGenomeGenerateRAM", settings["mem_mb"] * 1024 * 1024,
             "--outFileNamePrefix", str(Path(output) / "build.")], log)


def trim_command(r1, r2, out1, out2, json_report, html_report, threads, settings):
    command = ["fastp", "--in1", r1, "--in2", r2, "--out1", out1, "--out2", out2,
               "--thread", min(threads, 16), "--cut_tail", "--cut_tail_window_size", settings["quality_window"],
               "--cut_tail_mean_quality", settings["quality_mean"], "--length_required", settings["min_length"],
               "--json", json_report, "--html", html_report]
    if settings["adapter_r1"]:
        command += ["--adapter_sequence", settings["adapter_r1"]]
    if settings["adapter_r2"]:
        command += ["--adapter_sequence_r2", settings["adapter_r2"]]
    if not settings["adapter_r1"] and not settings["adapter_r2"]:
        command += ["--detect_adapter_for_pe"]
    if settings["poly_g"] == "yes":
        command += ["--trim_poly_g"]
    elif settings["poly_g"] == "no":
        command += ["--disable_trim_poly_g"]
    # No base correction, deduplication, merging, or low-complexity filtering.
    return command


def trim(r1, r2, out1, out2, json_report, html_report, threads, settings, logfile):
    parents(out1, out2, json_report, html_report, logfile)
    with open(logfile, "wb") as log:
        run(trim_command(r1, r2, out1, out2, json_report, html_report, threads, settings), log)


def check_index(index, fasta):
    lengths = {}
    with open(str(fasta) + ".fai") as stream:
        for line in stream:
            fields = line.split("\t")
            lengths[fields[0]] = int(fields[1])
    index = Path(index)
    names = (index / "chrName.txt").read_text().splitlines()
    sizes = (index / "chrLength.txt").read_text().splitlines()
    if len(names) != len(sizes) or dict(zip(names, map(int, sizes))) != lengths:
        raise ValueError("STAR index contigs/lengths do not match genome_fasta; rebuild the index")


def star_command(r1, r2, index, sample, prefix, threads, settings):
    groups = []
    for path in r1:
        run_id = Path(path).name.removesuffix("_1.fastq.gz")
        if groups:
            groups.append(",")
        groups += [f"ID:{run_id}", f"SM:{sample}", f"LB:{sample}", f"PU:{run_id}"]
    return ["STAR", "--runThreadN", max(1, threads - 1), "--genomeDir", index,
            "--readFilesIn", ",".join(r1), ",".join(r2), "--readFilesCommand", "zcat",
            "--outFileNamePrefix", prefix, "--outStd", "BAM_Unsorted", "--outSAMtype", "BAM", "Unsorted",
            "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD",
            "--outSAMattrRGline", *groups,
            "--outSAMprimaryFlag", "OneBestScore", "--outSAMmultNmax", "1",
            "--outSAMunmapped", "None",
            "--alignEndsType", "Local", "--outFilterMismatchNmax", "999",
            "--outFilterMismatchNoverLmax", "1",  # override STAR's second mismatch-fraction gate
            "--outFilterMismatchNoverReadLmax", settings["max_mismatch_fraction"],
            "--outFilterScoreMinOverLread", settings["min_score_fraction"],
            "--outFilterMatchNminOverLread", settings["min_match_fraction"],
            "--outFilterMultimapNmax", settings["multimap_max"],
            "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000"]


def align(r1, r2, index, fasta, sample, output, star_log, threads, settings, logfile):
    if threads < 2:
        raise ValueError("Alignment requires at least 2 allocated threads")
    check_index(index, fasta)
    parents(output, star_log, logfile)
    scratch = settings.get("scratch")
    if scratch:
        Path(scratch).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"insight-{sample}-", dir=scratch) as temporary:
        prefix = str(Path(temporary) / "star.")
        with open(logfile, "wb") as log:
            pipeline([star_command(r1, r2, index, sample, prefix, threads, settings),
                      ["samtools", "sort", "-@", "0", "-m", f'{settings["sort_mem_mb"]}M',
                       "-T", str(Path(temporary) / "sort"), "-o", output, "-"]], log)
            run(["samtools", "quickcheck", "-v", output], log)
        shutil.copyfile(prefix + "Log.final.out", star_log)


def finalize(source, fasta, sample, output, flagstat, idxstats, provenance, threads, settings, runs, logfile):
    parents(output, flagstat, idxstats, provenance, logfile)
    partial = str(output) + ".partial"
    try:
        with open(logfile, "wb") as log:
            with open(partial, "wb") as stream:
                pipeline([["samtools", "view", "-u", "-F", "0xB04", "-q", settings["min_mapq"], source],
                          ["samtools", "calmd", "-b", "-", fasta]], log, stdout=stream)
            run(["samtools", "quickcheck", "-v", partial], log)
            run(["samtools", "index", "-@", max(0, threads - 1), partial, str(output) + ".bai"], log)
            with open(flagstat, "wb") as stream:
                run(["samtools", "flagstat", partial], log, stdout=stream)
            # idxstats locates an index by BAM filename; use the explicit final path after promotion.
            Path(partial).replace(output)
            with open(idxstats, "wb") as stream:
                run(["samtools", "idxstats", output], log, stdout=stream)
            versions = {}
            for tool, args in (("STAR", ["--version"]), ("fastp", ["--version"]), ("samtools", ["--version"])):
                result = subprocess.run([tool] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
                versions[tool] = result.stdout.decode(errors="replace").splitlines()[0]
        Path(provenance).write_text(json.dumps({"sample": sample, "runs": list(runs), "settings": settings,
                                               "versions": versions, "reference": str(fasta),
                                               "coverage_filter": None, "duplicate_removal": False,
                                               "output": str(output)}, indent=2) + "\n")
    finally:
        Path(partial).unlink(missing_ok=True)
