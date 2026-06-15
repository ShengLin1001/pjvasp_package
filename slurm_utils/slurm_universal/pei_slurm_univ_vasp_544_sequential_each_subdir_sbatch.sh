#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
set -uo pipefail

# Walk the subdirectories below root_dir and submit each unfinished calculation
# as its own child Slurm job (sbatch --wait), one after another. This parent job
# only orchestrates; it does not run VASP itself, so it needs only 1 core and
# no VASP module (the child pei_slurm_univ_vasp_544.sh loads it).
# Workflow: this script -> pei_slurm_univ_submit --mode each-subdir -> sbatch --wait pei_slurm_univ_vasp_544.sh
#
# Usage (as a Slurm script):
#   sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh [root_dir] [lsubdir]
#   sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh "0.98,1.00,1.02"

_lib="$(command -v pei_slurm_univ_lib.sh || true)"
[[ -n "$_lib" ]] || { echo "❌ ERROR: pei_slurm_univ_lib.sh not found in PATH" >&2; exit 1; }
# shellcheck source=/dev/null
source "$_lib"

[[ -n "${SLURM_JOB_ID:-}" ]] || ps_fail "This is a Slurm batch script. Run it with: sbatch $(basename "$0")"

ROOT_DIR="./y_dir"
if [[ $# -ge 2 ]]; then ROOT_DIR="$1"; shift; fi
LSUBDIR="${1:-}"

echo "===== VASP 5.4.4 sequential (each-subdir sbatch) preflight ====="
ps_print_slurm_preflight
echo "📁 Sequential root: $ROOT_DIR    🧩 lsubdir: ${LSUBDIR:-all}    📄 child: pei_slurm_univ_vasp_544.sh"
echo "✅ Preflight passed. Starting sequential VASP (each-subdir sbatch)."
echo "================================================================"

exec pei_slurm_univ_submit --mode each-subdir \
    --root-dir "$ROOT_DIR" --lsubdir "$LSUBDIR" \
    --submit-script pei_slurm_univ_vasp_544.sh \
    --skip-if-file-contains OUTCAR "reached required accuracy"
