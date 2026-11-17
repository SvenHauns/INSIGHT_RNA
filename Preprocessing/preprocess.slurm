#!/usr/bin/env bash
#SBATCH --job-name=insight-preprocess
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm-preprocess-%j.log
#SBATCH --error=slurm-preprocess-%j.log

# Submit from INSIGHT_RNA/: sbatch Preprocessing/preprocess.slurm
# Supply your cluster's partition/account with sbatch options if required.
# This job runs FASTQ -> BAM preprocessing only; it never invokes main.py.
set -euo pipefail

# Edit these defaults, or export the variables before submitting.
project_dir=${INSIGHT_RNA_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}
config_file=${INSIGHT_PREPROCESS_CONFIG:-Preprocessing/config.json}
conda_env=${INSIGHT_CONDA_ENV-insight-preprocess}
use_conda=${INSIGHT_USE_CONDA:-1}

# If required by your cluster, load its Conda module here before activation.
# Alternatively set INSIGHT_CONDA_SH=/absolute/path/to/etc/profile.d/conda.sh.
if [[ -n "$conda_env" ]]; then
    if [[ -n ${INSIGHT_CONDA_SH:-} ]]; then
        source "$INSIGHT_CONDA_SH"
    elif command -v conda >/dev/null 2>&1; then
        conda_base=$(conda info --base)
        source "$conda_base/etc/profile.d/conda.sh"
    else
        printf '%s\n' 'Conda not found. Load its module or set INSIGHT_CONDA_SH before submitting.' >&2
        exit 1
    fi
    conda activate "$conda_env"
fi

cd -- "$project_dir"
if [[ ! -f preprocess.py ]]; then
    printf '%s\n' 'Submit from INSIGHT_RNA/, or export INSIGHT_RNA_DIR to its absolute path.' >&2
    exit 1
fi
if [[ ! -f "$config_file" ]]; then
    printf 'Missing config: %s. Copy/edit Preprocessing/config.example.json first.\n' "$config_file" >&2
    exit 1
fi
for tool in python3 snakemake; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        printf 'Required executable is not on PATH: %s\n' "$tool" >&2
        exit 1
    fi
done

cores=${SLURM_CPUS_PER_TASK:-8}
node_mem_mb=${SLURM_MEM_PER_NODE:-81920}
# Reserve memory for Snakemake and other process overhead outside rule requests.
workflow_mem_mb=$((node_mem_mb - 4096))
if (( cores < 2 || workflow_mem_mb < 1 )); then
    printf '%s\n' 'Allocate at least two CPUs and more than 4 GiB of node memory.' >&2
    exit 1
fi

# This template uses a JSON configuration, matching config.example.json.
python3 - "$config_file" "$cores" "$workflow_mem_mb" <<'PY'
import json
from pathlib import Path
import sys

if sys.version_info < (3, 10):
    raise SystemExit("Python 3.10+ is required")
settings = json.loads(Path(sys.argv[1]).read_text())
threads, memory = settings.get("threads", 8), settings.get("mem_mb", 64000)
if threads > int(sys.argv[2]):
    raise SystemExit("Config threads exceeds Slurm CPUs: increase --cpus-per-task or lower config threads")
if memory > int(sys.argv[3]):
    raise SystemExit("Config mem_mb exceeds Slurm allocation minus headroom: increase --mem or lower config mem_mb")
PY

extra_args=()
case "$use_conda" in
    1)
        if ! command -v conda >/dev/null 2>&1; then
            printf '%s\n' 'Conda must be available for rule environments (--use-conda).' >&2
            exit 1
        fi
        extra_args+=(--use-conda --conda-frontend "${INSIGHT_CONDA_FRONTEND:-conda}")
        ;;
    0)
        for tool in fastp STAR samtools; do
            if ! command -v "$tool" >/dev/null 2>&1; then
                printf 'Required preprocessing tool is not on PATH: %s\n' "$tool" >&2
                exit 1
            fi
        done
        ;;
    *) printf '%s\n' 'INSIGHT_USE_CONDA must be 0 or 1.' >&2; exit 2 ;;
esac
for option in "$@"; do
    case "$option" in
        --dry-run|--keep-temp) extra_args+=("$option") ;;
        *) printf 'Unsupported job argument: %s (allowed: --dry-run, --keep-temp)\n' "$option" >&2; exit 2 ;;
    esac
done

printf 'Job: %s | Host: %s | CPUs: %s | Workflow memory: %s MB\n' \
    "${SLURM_JOB_ID:-local}" "$(hostname)" "$cores" "$workflow_mem_mb"
printf 'Project: %s\nConfig: %s\nStarted: %s\n' "$PWD" "$config_file" "$(date -Is)"

exec python3 -u preprocess.py --config "$config_file" \
    --cores "$cores" --mem-mb "$workflow_mem_mb" \
    --scheduler "${INSIGHT_SCHEDULER:-greedy}" "${extra_args[@]}"
