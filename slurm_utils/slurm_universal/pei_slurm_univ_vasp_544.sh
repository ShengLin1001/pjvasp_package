#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out
set -uo pipefail

# Run or resume one VASP 5.4.4 calculation in the submission directory.
# Workflow: pei_slurm_univ_vasp_544.sh -> pei_vasp_univ_sbatch

_lib="$(command -v pei_slurm_univ_lib.sh || true)"
[[ -n "$_lib" ]] || { echo "❌ ERROR: pei_slurm_univ_lib.sh not found in PATH" >&2; exit 1; }
_env="$(command -v pei_vasp_univ_load_env || true)"
[[ -n "$_env" ]] || { echo "❌ ERROR: pei_vasp_univ_load_env not found in PATH" >&2; exit 1; }
# shellcheck source=/dev/null
source "$_lib"
# shellcheck source=/dev/null
source "$_env"

[[ -n "${SLURM_JOB_ID:-}" ]] || ps_fail "This is a Slurm batch script. Run it with: sbatch $(basename "$0")"

echo "================ VASP 5.4.4 Slurm preflight ================"
ps_print_slurm_preflight
pei_vasp_univ_load_env
pei_vasp_check_inputs INCAR POSCAR POTCAR KPOINTS Y_CONSTR_LATT
echo "✅ Preflight passed. Running single-directory VASP workflow."
echo "==========================================================="

pei_vasp_univ_sbatch "." "$VASP_EXE"
status=$?
(( status == 10 )) && { ps_ok "Already completed; nothing to do."; exit 0; }
exit "$status"
