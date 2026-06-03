#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

MODULE_FILE="/public3/soft/modules/module.sh"
VASP_BIN_DIR="/public3/home/scg6928/mysoft/vasp/vasp/544-yin/vasp.5.4.4.pl2/bin"
VASP_EXE="vasp_std"
ROOT_DIR="./y_dir"
SEQUENTIAL_HELPER="pei_vasp_univ_sbatch_sequential"

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

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    fail "This is a Slurm batch script. Please run it with: sbatch sub.544.sequential.sh"
fi

LSUBDIR="${1:-}"

echo "============ VASP 5.4.4 sequential Slurm preflight ==========="
echo "📄 Slurm output: slurm-${SLURM_JOB_ID}.out"
echo "📁 Work dir: $(pwd)"
echo "📁 Sequential root: $ROOT_DIR"
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
command -v "$VASP_EXE" >/dev/null 2>&1 || fail "VASP executable not found in PATH: $VASP_EXE"
ok "VASP executable: $(command -v "$VASP_EXE")"

[[ -d "$ROOT_DIR" ]] || fail "sequential root directory not found in $(pwd): $ROOT_DIR"
ok "Found sequential root directory: $ROOT_DIR"

SEQUENTIAL_HELPER_PATH="$(command -v "$SEQUENTIAL_HELPER" || true)"
[[ -n "$SEQUENTIAL_HELPER_PATH" ]] || fail "sequential helper not found in PATH: $SEQUENTIAL_HELPER"
[[ -x "$SEQUENTIAL_HELPER_PATH" ]] || fail "sequential helper is not executable: $SEQUENTIAL_HELPER_PATH"
ok "Sequential helper: $SEQUENTIAL_HELPER_PATH"

if [[ $# -gt 0 ]]; then
    warn "Extra arguments passed to sequential helper after root and executable: $*"
fi

echo "✅ Preflight passed. Starting sequential VASP workflow."
echo "============================================================="
"$SEQUENTIAL_HELPER_PATH" "$ROOT_DIR" "$VASP_EXE" "$@"
