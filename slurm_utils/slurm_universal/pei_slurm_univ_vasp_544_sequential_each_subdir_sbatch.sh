#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

MODULE_FILE="/public3/soft/modules/module.sh"
VASP_BIN_DIR="/public3/home/scg6928/mysoft/vasp/vasp/544-yin/vasp.5.4.4.pl2/bin"
ROOT_DIR="./y_dir"
SUBMIT_SCRIPT="pei_slurm_univ_vasp_544.sh"
SEQUENTIAL_HELPER="pei_vasp_univ_sbatch_sequential_each_subdir_sbatch"

ok() {
    echo "✅ $*"
}

warn() {
    echo "⚠️  $*"
}

fail() {
    echo "❌ ERROR: $*" >&2
    exit 1
}

check_parent_resources() {
    local nodes="${SLURM_JOB_NUM_NODES:-}"
    local tasks="${SLURM_NTASKS:-}"

    [[ "$nodes" == "1" ]] || fail "sequential parent job must request exactly 1 node; got SLURM_JOB_NUM_NODES=${nodes:-unset}"
    [[ "$tasks" == "1" ]] || fail "sequential parent job must request exactly 1 task/core; got SLURM_NTASKS=${tasks:-unset}"
}

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    fail "This is a Slurm batch script. Please run it with: sbatch pei_slurm_univ_vasp_544_sequential_each_subdir_sbatch.sh"
fi
check_parent_resources

if [[ $# -ge 2 ]]; then
    ROOT_DIR="$1"
    shift
fi
LSUBDIR="${1:-}"

echo "============ VASP 5.4.4 sequential Slurm preflight ==========="
echo "📄 Slurm output: slurm-${SLURM_JOB_ID}.out"
echo "📁 Work dir: $(pwd)"
echo "📁 Sequential root: $ROOT_DIR"
echo "📄 Child submit script: $SUBMIT_SCRIPT"
echo "🧩 lsubdir: ${LSUBDIR:-all}"
echo "🖥️  Host: $(hostname)"
echo "🧾 Job ID: ${SLURM_JOB_ID}"
echo "📦 Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo "🔢 Nodes: ${SLURM_JOB_NUM_NODES:-unknown}"
echo "🔢 Tasks: ${SLURM_NTASKS:-unknown}"

[[ -r "$MODULE_FILE" ]] || fail "module file not readable: $MODULE_FILE"
source "$MODULE_FILE" || fail "failed to source module file: $MODULE_FILE"
module load mpi/intel/17.0.7-thc || fail "failed to load module: mpi/intel/17.0.7-thc"
ok "Loaded module: mpi/intel/17.0.7-thc"

[[ -d "$VASP_BIN_DIR" ]] || fail "VASP bin directory not found: $VASP_BIN_DIR"
export PATH="$VASP_BIN_DIR:$PATH"
command -v sbatch >/dev/null 2>&1 || fail "sbatch command not found in PATH"
ok "sbatch executable: $(command -v sbatch)"

[[ -d "$ROOT_DIR" ]] || fail "sequential root directory not found in $(pwd): $ROOT_DIR"
ok "Found sequential root directory: $ROOT_DIR"

SEQUENTIAL_HELPER_PATH="$(command -v "$SEQUENTIAL_HELPER" || true)"
[[ -n "$SEQUENTIAL_HELPER_PATH" ]] || fail "sequential helper not found in PATH: $SEQUENTIAL_HELPER"
[[ -x "$SEQUENTIAL_HELPER_PATH" ]] || fail "sequential helper is not executable: $SEQUENTIAL_HELPER_PATH"
ok "Sequential helper: $SEQUENTIAL_HELPER_PATH"

if [[ $# -gt 1 ]]; then
    warn "Extra arguments passed to sequential helper after lsubdir: ${*:2}"
fi

echo "✅ Preflight passed. Starting sequential VASP workflow."
echo "============================================================="
"$SEQUENTIAL_HELPER_PATH" "$ROOT_DIR" "$SUBMIT_SCRIPT" "$@"
