#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

GLOBAL = {
    "samples": "data/samples.tsv",
    "oops_samples": "data/samples.tsv",
    "nc": "data/samples.tsv",
    "cl": "data/samples.tsv",
    "reference": "ref/genome.fa",
    "mapped": "Preprocessing/mapped/",
    "gff_path": "refseq_with_start_codon",
    "genome": "hg38",
    "run_type": "3UTR",
    "deseq_output": "deseq_output/",
    "exons": "files/exons.bed",
    "model_path": "models/model_500_2_small.pt",
    "eclip": "eClip/",
}

PIPELINES = [
    {
        "name": "OOPS-seq",
        "snakefile": "OOPS/Snakefile",
        "extra_config": {
            "filter_threshold": 0.05,
        },
    },
    {
        "name": "DMS",
        "snakefile": "DMS/Snakefile",
        "extra_config": {
            "bit_data": "hg38.2bit"
        },
    },
    {
        "name": "Integration",
        "snakefile": "Integration/Snakefile",
        "extra_config": {
        },
    },
]



def as_config_args(d: dict):
    items = []
    for k, v in d.items():
        items.append(f"{k}={v}")
    return ["--config"] + items if items else []

def run_pipeline(snakefile: str, name: str, config: dict, cores: int, dry_run: bool, use_conda: bool, profile: str | None):

    cmd = [
        "snakemake",
        "-s", snakefile,
        "--cores", str(cores),
        "--rerun-incomplete",
        "--printshellcmds",
        "--nolock",
    ]
    if dry_run:
        cmd.append("--dry-run")
    if use_conda:
        cmd.append("--use-conda")
    if profile:
        cmd += ["--profile", profile]

    cmd += as_config_args(config)

    print(f"\n=== Running {name} ===")
    print(" ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
        print(f"=== {name}: OK ===\n")
    except subprocess.CalledProcessError as e:
        print(f"*** {name}: FAILED with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)

def main():
    parser = argparse.ArgumentParser(description="Run three snakemake pipelines sequentially with shared Globals.")
    parser.add_argument("--cores", type=int, default=1, help="cores to pass to each pipeline.")
    parser.add_argument("--dry-run", action="store_true", help="snakemake --dry-run for all pipelines.")
    parser.add_argument("--use-conda", action="store_true", help="enable snakemake --use-conda.")
    parser.add_argument("--profile", type=str, default=None, help="Snakemake --profile name/path (optional).")
    parser.add_argument("--samples_oops", type=str, default=None, help="samples to process")
    parser.add_argument("--nc_samples_oops", type=str, default=None, help="samples to process")
    parser.add_argument("--cl_samples_oops", type=str, default=None, help="samples to process")
    parser.add_argument("--samples", type=str, default=None, help="samples to process")
    parser.add_argument("--run_type", type=str, default=None, help="DMS run type, 3UTR, 5UTR")
    parser.add_argument("--model_path", type=str, default=None, help="path to RibonanzaNet model")
    
    args = parser.parse_args()
    
    if args.samples_oops != None: GLOBAL["oops_samples"] = args.samples_oops
    if args.samples_oops != None: GLOBAL["nc"] = args.nc_samples_oops
    if args.samples_oops != None: GLOBAL["cl"] = args.cl_samples_oops
    if args.samples != None: GLOBAL["samples"] = args.samples
    if args.run_type != None: GLOBAL["run_type"] = args.run_type
    if args.model_path != None: GLOBAL["model_path"] = args.model_path
    
    
    for p in PIPELINES:
        merged_cfg = {**GLOBAL, **p.get("extra_config", {})}
        run_pipeline(
            snakefile=p["snakefile"],
            name=p["name"],
            config=merged_cfg,
            cores=args.cores,
            dry_run=args.dry_run,
            use_conda=args.use_conda,
            profile=args.profile,
        )

    print("All pipelines finished successfully.")

if __name__ == "__main__":
    main()
