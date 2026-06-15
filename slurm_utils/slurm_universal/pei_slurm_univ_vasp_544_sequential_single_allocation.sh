#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
set -uo pipefail

# Run unfinished VASP calculations below root_dir one at a time, inside this
# single Slurm allocation (each calculation runs via srun here).
# Workflow: this script -> pei_slurm_univ_submit --mode single-alloc -> pei_vasp_univ_sbatch
#
# Usage (as a Slurm script):
#   sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh [root_dir] [lsubdir]
#   sbatch pei_slurm_univ_vasp_544_sequential_single_allocation.sh "0.98,1.00,1.02"

_lib="$(command -v pei_slurm_univ_lib.sh || true)"
[[ -n "$_lib" ]] || { echo "❌ ERROR: pei_slurm_univ_lib.sh not found in PATH" >&2; exit 1; }
_env="$(command -v pei_vasp_univ_load_env || true)"
[[ -n "$_env" ]] || { echo "❌ ERROR: pei_vasp_univ_load_env not found in PATH" >&2; exit 1; }
# shellcheck source=/dev/null
source "$_lib"
# shellcheck source=/dev/null
source "$_env"

[[ -n "${SLURM_JOB_ID:-}" ]] || ps_fail "This is a Slurm batch script. Run it with: sbatch $(basename "$0")"

ROOT_DIR="./y_dir"
if [[ $# -ge 2 ]]; then ROOT_DIR="$1"; shift; fi
LSUBDIR="${1:-}"

echo "===== VASP 5.4.4 sequential (single allocation) preflight ====="
ps_print_slurm_preflight
echo "📁 Sequential root: $ROOT_DIR    🧩 lsubdir: ${LSUBDIR:-all}"
pei_vasp_univ_load_env
echo "✅ Preflight passed. Starting sequential VASP (single allocation)."
echo "==============================================================="

exec pei_slurm_univ_submit --mode single-alloc \
    --root-dir "$ROOT_DIR" --lsubdir "$LSUBDIR" \
    --run-cmd "pei_vasp_univ_sbatch . $VASP_EXE"
