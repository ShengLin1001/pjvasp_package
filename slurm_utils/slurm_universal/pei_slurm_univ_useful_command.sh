#!/bin/bash

if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
    reset="$(printf '\033[0m')"
    bold="$(printf '\033[1m')"
    title_color="$(printf '\033[1;36m')"
    section_color="$(printf '\033[1;35m')"
    command_color="$(printf '\033[0;32m')"
    desc_color="$(printf '\033[0;37m')"
else
    reset=""
    bold=""
    title_color=""
    section_color=""
    command_color=""
    desc_color=""
fi

section() {
    printf '\n%s[%s]%s\n' "$section_color" "$1" "$reset"
}

command_item() {
    printf '  %s%s%s\n' "$command_color" "$1" "$reset"
    printf '      %s%s%s\n\n' "$desc_color" "$2" "$reset"
}

printf '%s%s%s\n' "$title_color" "Useful commands" "$reset"
printf '%s%s%s\n' "$bold" "===============" "$reset"

section "Slurm"
command_item "sinfo" \
    "Show partition and node status."
command_item 'squeue -u "$USER"' \
    "Show jobs submitted by the current user."
command_item "scontrol show job <jobid>" \
    "Show detailed information for a job."
command_item "scancel <jobid>" \
    "Cancel a job."
command_item "sacct -j <jobid> --format=JobID,JobName,Partition,State,Elapsed,MaxRSS,ExitCode" \
    "Show accounting information for a finished job."
command_item "tail -f slurm-<jobid>.out" \
    "Follow the Slurm output file."

section "VASP"
command_item "head -n 10 POTCAR" \
    "Show the first lines of POTCAR."
command_item "cat -n POSCAR" \
    "Print POSCAR with line numbers."
command_item 'grep -A 5 "TOTAL-FORCE" OUTCAR' \
    "Show force blocks in OUTCAR. Use atom_count + 3 lines when needed."
command_item 'grep "in kB" OUTCAR' \
    "Show stress tensor lines in OUTCAR."

section "VS Code"
command_item "Ctrl + Shift + P" \
    "Open the command palette."
command_item "Ctrl + /" \
    "Comment or uncomment the current line or selection."

section "Codex"
command_item "Ctrl + N" \
    "Start a new chat."

section "Network"
command_item "curl -I https://www.baidu.com" \
    "Check whether the site is reachable and print response headers."
command_item "curl ipinfo.io" \
    "Show public IP information."
command_item "curl cip.cc" \
    "Show public IP and location information."
