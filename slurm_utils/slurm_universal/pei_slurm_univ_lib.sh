#!/bin/bash
# Shared helpers for the pei_slurm_univ_* batch submission scripts.
#
# Source this file; do not execute it:
#
#   source "$(command -v pei_slurm_univ_lib.sh)"
#
# Provides:
#   ps_log / ps_ok / ps_warn / ps_fail   informative output (ps_fail exits 1)
#   ps_print_slurm_preflight             echo the Slurm job-info block
#   ps_resolve_job_dirs <root> <lsubdir> fill global array PS_JOB_DIRS
#   ps_resolve_script_path <script>      print absolute path (file or PATH lookup)
#   ps_report_failures <header> [dir...] print a header + the listed dirs

# --- informative output ------------------------------------------------------
ps_log()  { echo "$*"; }
ps_ok()   { echo "✅ $*"; }
ps_warn() { echo "⚠️  $*"; }
ps_fail() { echo "❌ ERROR: $*" >&2; exit 1; }

# --- Slurm preflight echo block ----------------------------------------------
ps_print_slurm_preflight() {
    echo "📄 Slurm output: slurm-${SLURM_JOB_ID:-NA}.out"
    echo "📁 Work dir: $(pwd)"
    echo "🖥️  Host: $(hostname)"
    echo "🧾 Job ID: ${SLURM_JOB_ID:-unknown}"
    echo "📦 Partition: ${SLURM_JOB_PARTITION:-unknown}"
    echo "🔢 Nodes: ${SLURM_JOB_NUM_NODES:-unknown}"
    echo "🔢 Tasks: ${SLURM_NTASKS:-unknown}"
}

# --- resolve the directories to process --------------------------------------
# Fills the global array PS_JOB_DIRS with first-level subdirectories of <root>.
# When <lsubdir> is a non-empty comma-separated basename list, only those are
# used; otherwise every first-level subdirectory is used.
ps_resolve_job_dirs() {
    local root_dir="$1" lsubdir="$2"
    PS_JOB_DIRS=()
    [[ -d "$root_dir" ]] || ps_fail "root directory not found: $root_dir"
    # nullglob so an empty directory yields an empty array, not a literal glob.
    shopt -s nullglob
    if [[ -z "$lsubdir" ]]; then
        PS_JOB_DIRS=("$root_dir"/*/)
    else
        local subdir job_dir
        local -a items
        IFS=',' read -r -a items <<< "$lsubdir"
        for subdir in "${items[@]}"; do
            [[ -n "$subdir" ]]     || ps_fail "empty directory name in lsubdir: $lsubdir"
            [[ "$subdir" != */* ]] || ps_fail "lsubdir entries must be basenames below $root_dir: $subdir"
            job_dir="${root_dir%/}/$subdir/"
            [[ -d "$job_dir" ]]    || ps_fail "lsubdir directory not found: $job_dir"
            PS_JOB_DIRS+=("$job_dir")
        done
    fi
    (( ${#PS_JOB_DIRS[@]} > 0 )) || ps_fail "no job directories found below: $root_dir"
}

# --- resolve a submit script to an absolute path -----------------------------
# Accepts either a file path (containing '/' or an existing file) or a bare
# command name to look up on PATH. Prints the absolute, readable path.
ps_resolve_script_path() {
    local script="$1" path
    if [[ "$script" == */* || -f "$script" ]]; then
        [[ -f "$script" ]] || ps_fail "submit script not found: $script"
        path="$(cd -- "$(dirname -- "$script")" && pwd)/$(basename -- "$script")"
    else
        path="$(command -v "$script" || true)"
    fi
    [[ -n "$path" ]] || ps_fail "submit script not found as file or command in PATH: $script"
    [[ -r "$path" ]] || ps_fail "submit script is not readable: $path"
    printf '%s\n' "$path"
}

# --- report a list of directories under a header -----------------------------
ps_report_failures() {
    local header="$1"; shift
    (( $# == 0 )) && return 0
    echo "$header"
    local d
    for d in "$@"; do
        echo "  - 📁 $d"
    done
}
