#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.out

MODULE_FILE="/public3/soft/modules/module.sh"
VASP_BIN_DIR="/public3/home/scg6928/mysoft/vasp/vasp/544-yin/vasp.5.4.4.pl2/bin"
VASP_EXE="vasp_std"
SINGLE_JOB_RUNNER="pei_vasp_univ_sbatch"

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
    fail "This is a Slurm batch script. Please run it with: sbatch sub.544.sh"
fi

echo "================ VASP 5.4.4 Slurm preflight ================"
echo "📄 Slurm output: slurm-${SLURM_JOB_ID}.out"
echo "📁 Work dir: $(pwd)"
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

SINGLE_JOB_RUNNER_PATH="$(command -v "$SINGLE_JOB_RUNNER" || true)"
[[ -n "$SINGLE_JOB_RUNNER_PATH" ]] || fail "single job runner not found in PATH: $SINGLE_JOB_RUNNER"
[[ -x "$SINGLE_JOB_RUNNER_PATH" ]] || fail "single job runner is not executable: $SINGLE_JOB_RUNNER_PATH"
ok "Single job runner: $SINGLE_JOB_RUNNER_PATH"

for input_file in INCAR POSCAR POTCAR KPOINTS Y_CONSTR_LATT; do
    [[ -f "$input_file" ]] || fail "required VASP input file missing in $(pwd): $input_file"
    ok "Found input file: $input_file"
done

# 或许后面可以自动生成POTCAR利用vaspkit

if [[ $# -gt 0 ]]; then
    warn "Extra arguments passed to single job runner after job directory and executable: $*"
fi

echo "✅ Preflight passed. Starting single-directory VASP workflow."
echo "============================================================="
"$SINGLE_JOB_RUNNER_PATH" "." "$VASP_EXE" "$@"
status=$?

if (( status == 10 )); then
    ok "Single-directory workflow skipped an already completed calculation."
    exit 0
fi

exit "$status"
