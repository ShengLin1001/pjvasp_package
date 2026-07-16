#!/bin/bash
#
# pei_slurm_univ_useful_command.sh
#
# Without options: print the useful-command cheat sheet (Slurm / VASP / editor / network).
# With options: actually probe network reachability of the endpoints that matter on this
# HPC login node (OpenAI, Claude, GitHub, PyPI, conda), show the public IP, or dump the
# proxy environment. Handy right after setting up a tunnel or a new proxy.
#
# Examples:
#   pei_slurm_univ_useful_command.sh                  # cheat sheet only
#   pei_slurm_univ_useful_command.sh --claude         # probe Anthropic / Claude endpoints
#   pei_slurm_univ_useful_command.sh --openai --proxy # probe OpenAI + show proxy vars
#   pei_slurm_univ_useful_command.sh --net --ip       # probe every group + public IP
#   pei_slurm_univ_useful_command.sh --monitor        # summarize running Slurm jobs
#   pei_slurm_univ_useful_command.sh --all --timeout 3

script_name="$(basename "$0")"
timeout_default=8
monitor_interval_default=5
stale_after_default=600

################ colors ################

if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
    reset="$(printf '\033[0m')"
    bold="$(printf '\033[1m')"
    title_color="$(printf '\033[1;36m')"
    section_color="$(printf '\033[1;35m')"
    command_color="$(printf '\033[0;32m')"
    desc_color="$(printf '\033[0;37m')"
    ok_color="$(printf '\033[0;32m')"
    warn_color="$(printf '\033[0;33m')"
    fail_color="$(printf '\033[0;31m')"
else
    reset=""
    bold=""
    title_color=""
    section_color=""
    command_color=""
    desc_color=""
    ok_color=""
    warn_color=""
    fail_color=""
fi
# to here

################ helpers ################

section() {
    printf '\n%s[%s]%s\n' "$section_color" "$1" "$reset"
}

command_item() {
    printf '  %s%s%s\n' "$command_color" "$1" "$reset"
    printf '      %s%s%s\n\n' "$desc_color" "$2" "$reset"
}

fail() {
    printf '%s❌ ERROR: %s%s\n' "$fail_color" "$1" "$reset" >&2
    exit 1
}

warn() {
    printf '  %s⚠️  %s%s\n' "$warn_color" "$1" "$reset"
}

# Structural check only: the timeout has to be a positive integer, curl rejects the rest.
check_positive_int() {
    local value="$1"
    local name="$2"
    case "$value" in
        ''|*[!0-9]*) fail "$name must be a positive integer, got '$value'" ;;
    esac
    [ "$value" -gt 0 ] || fail "$name must be a positive integer, got '$value'"
    printf '%s\n' "$value"
}
# to here

################ endpoint table ################

# Keys are the CLI flags (--openai, --claude, ...); values are space separated URLs.
# Endpoints are picked so that a single group answers one practical question:
# "can I use this tool from this node?"
declare -A dict_group_url=(
    [base]="https://www.baidu.com/ https://www.google.com/ https://1.1.1.1/"
    [openai]="https://api.openai.com/ https://api.openai.com/v1/models https://chatgpt.com/"
    [claude]="https://api.anthropic.com/ https://api.anthropic.com/v1/messages https://claude.ai/"
    [github]="https://github.com/ https://api.github.com/ https://raw.githubusercontent.com/"
    [pypi]="https://pypi.org/simple/ https://files.pythonhosted.org/ https://pypi.tuna.tsinghua.edu.cn/simple/"
    [conda]="https://repo.anaconda.com/pkgs/main/ https://conda.anaconda.org/conda-forge/ https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/"
)

declare -A dict_group_desc=(
    [base]="Baseline reachability: domestic, overseas, plain IP (DNS vs routing)."
    [openai]="OpenAI API and web app (Codex / ChatGPT)."
    [claude]="Anthropic API and web app (Claude Code / claude.ai)."
    [github]="GitHub web, API and raw download (git clone, pip install from git)."
    [pypi]="PyPI index, package CDN and the Tsinghua mirror (pip install)."
    [conda]="Anaconda / conda-forge channels and the Tsinghua mirror (conda install)."
)

# Assoc arrays have no stable order, so keep the display order explicit.
lgroup_order=(base openai claude github pypi conda)
# to here

################ network probes ################

n_ok=0
n_warn=0
n_fail=0
lfailed=()

# HEAD is cheap, but several endpoints (raw.githubusercontent.com, some mirrors, API
# routes that only accept POST) answer it with 405/403/000 even though they are perfectly
# reachable. Retry those with a 1-byte ranged GET so a picky verb is not reported as
# "network down". Prints: "<http_code> <seconds> <method>".
net_probe_url() {
    local url="$1"
    local out code elapsed method

    method="HEAD"
    out="$(curl -s -o /dev/null -I --max-time "$timeout" \
        -w '%{http_code} %{time_total}' "$url" 2>/dev/null)"
    code="${out%% *}"
    elapsed="${out##* }"

    if [ -z "$code" ] || [ "$code" = "000" ] || [ "$code" = "405" ] || [ "$code" = "501" ]; then
        method="GET"
        out="$(curl -s -o /dev/null -r 0-0 --max-time "$timeout" \
            -w '%{http_code} %{time_total}' "$url" 2>/dev/null)"
        code="${out%% *}"
        elapsed="${out##* }"
    fi

    [ -z "$code" ] && code="000"
    [ -z "$elapsed" ] && elapsed="0.000"
    printf '%s %s %s\n' "$code" "$elapsed" "$method"
}

# What we test is the network path, not the API contract: 401/403/404/421 all prove the
# request reached the server (we deliberately send no API key), so any HTTP code counts as
# ✅ reachable. Only 000 (DNS / TCP / TLS / timeout: no answer at all) is a real ❌, and
# 5xx gets a ⚠️ because there the path works but the far end is unhealthy.
net_check_group() {
    local group="$1"
    local url result code elapsed method

    section "Network / --${group}"
    printf '  %s%s%s\n\n' "$desc_color" "${dict_group_desc[$group]}" "$reset"

    for url in ${dict_group_url[$group]}; do
        result="$(net_probe_url "$url")"
        read -r code elapsed method <<< "$result"

        if [ "$code" = "000" ]; then
            printf '  %s❌ %-3s%s  %-52s  %ss  %s\n' \
                "$fail_color" "$code" "$reset" "$url" "$elapsed" "$method"
            lfailed+=("$url")
            n_fail=$((n_fail + 1))
        elif [ "$code" -ge 500 ] 2>/dev/null; then
            printf '  %s⚠️  %-3s%s  %-52s  %ss  %s\n' \
                "$warn_color" "$code" "$reset" "$url" "$elapsed" "$method"
            n_warn=$((n_warn + 1))
        else
            printf '  %s✅ %-3s%s  %-52s  %ss  %s\n' \
                "$ok_color" "$code" "$reset" "$url" "$elapsed" "$method"
            n_ok=$((n_ok + 1))
        fi
    done
}

net_show_ip() {
    local service out

    section "Network / --ip"
    for service in "ipinfo.io" "cip.cc"; do
        printf '  %s▶️  curl %s%s\n' "$command_color" "$service" "$reset"
        out="$(curl -s --max-time "$timeout" "$service" 2>/dev/null)"
        if [ -z "$out" ]; then
            warn "no answer from $service within ${timeout}s"
        else
            printf '%s\n' "$out" | while IFS= read -r line; do
                printf '      %s%s%s\n' "$desc_color" "$line" "$reset"
            done
        fi
        printf '\n'
    done
}

net_show_proxy() {
    local lvar=(http_proxy https_proxy all_proxy no_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY NO_PROXY)
    local var value found=0

    section "Network / --proxy"
    for var in "${lvar[@]}"; do
        value="${!var:-}"
        if [ -n "$value" ]; then
            printf '  %s✅ %-12s%s = %s\n' "$ok_color" "$var" "$reset" "$value"
            found=1
        fi
    done
    if [ "$found" -eq 0 ]; then
        warn "no proxy variable is set; traffic goes out directly"
        printf '      %sexport https_proxy=http://127.0.0.1:7890 http_proxy=http://127.0.0.1:7890%s\n' \
            "$desc_color" "$reset"
    fi
}

print_summary() {
    local url

    printf '\n%s================ 📊 summary%s\n' "$bold" "$reset"
    printf '  %s✅ reachable: %d%s   %s⚠️  server error (5xx): %d%s   %s❌ unreachable: %d%s\n' \
        "$ok_color" "$n_ok" "$reset" \
        "$warn_color" "$n_warn" "$reset" \
        "$fail_color" "$n_fail" "$reset"
    printf '  %s401/403/404/421 still count as reachable: no API key is sent, only the path is tested.%s\n' \
        "$desc_color" "$reset"

    if [ "$n_fail" -gt 0 ]; then
        printf '\n  %sUnreachable endpoints:%s\n' "$fail_color" "$reset"
        for url in "${lfailed[@]}"; do
            printf '    %s❌ %s%s\n' "$fail_color" "$url" "$reset"
        done
        return 1
    fi

    printf '\n  %s🎉 all probed endpoints answered%s\n' "$ok_color" "$reset"
    return 0
}
# to here

################ Slurm monitor ################

get_job_field() {
    local job_info="$1"
    local field="$2"
    local token

    for token in $job_info; do
        case "$token" in
            "$field"=*)
                printf '%s\n' "${token#*=}"
                return 0
                ;;
        esac
    done
    return 1
}

resolve_job_path() {
    local path_raw="$1"
    local path_workdir="$2"

    case "$path_raw" in
        ''|'(null)'|'N/A'|'Unknown')
            printf '\n'
            ;;
        /*)
            printf '%s\n' "$path_raw"
            ;;
        *)
            printf '%s/%s\n' "${path_workdir%/}" "$path_raw"
            ;;
    esac
}

is_monitorable_path() {
    local path_file="$1"
    [ -n "$path_file" ] && [ "$path_file" != "/dev/null" ]
}

snapshot_file() {
    local path_file="$1"
    local snapshot

    if is_monitorable_path "$path_file" && [ -e "$path_file" ]; then
        snapshot="$(stat -c '%s %Y' -- "$path_file" 2>/dev/null)"
        if [ -n "$snapshot" ]; then
            printf '1 %s\n' "$snapshot"
            return 0
        fi
    fi
    printf '0 0 0\n'
}

duration_to_seconds() {
    local duration="$1"
    local days=0
    local clock="$duration"
    local first=0
    local second=0
    local third=0
    local n_seconds=0
    local lpart=()

    if [[ "$clock" == *-* ]]; then
        days="${clock%%-*}"
        clock="${clock#*-}"
    fi
    IFS=':' read -r -a lpart <<< "$clock"
    case "${#lpart[@]}" in
        3)
            first="${lpart[0]}"
            second="${lpart[1]}"
            third="${lpart[2]}"
            ;;
        2)
            second="${lpart[0]}"
            third="${lpart[1]}"
            ;;
        *)
            printf '0\n'
            return 0
            ;;
    esac

    case "$days:$first:$second:$third" in
        *[!0-9:]*) printf '0\n'; return 0 ;;
    esac
    n_seconds=$((10#$days * 86400 + 10#$first * 3600 + 10#$second * 60 + 10#$third))
    printf '%d\n' "$n_seconds"
}

format_seconds() {
    local n_seconds="$1"
    local days hours minutes seconds

    [ "$n_seconds" -ge 0 ] 2>/dev/null || n_seconds=0
    days=$((n_seconds / 86400))
    hours=$(((n_seconds % 86400) / 3600))
    minutes=$(((n_seconds % 3600) / 60))
    seconds=$((n_seconds % 60))
    if [ "$days" -gt 0 ]; then
        printf '%dd %02d:%02d:%02d' "$days" "$hours" "$minutes" "$seconds"
    else
        printf '%02d:%02d:%02d' "$hours" "$minutes" "$seconds"
    fi
}

# A short quiet sample is normal for buffered programs. Only mark it stale after the
# configured age threshold, and keep the wording advisory rather than definitive.
classify_file_activity() {
    local before_exists="$1"
    local before_size="$2"
    local before_mtime="$3"
    local after_exists="$4"
    local after_size="$5"
    local after_mtime="$6"
    local now="$7"
    local stale_after="$8"

    file_state="MISSING"
    file_delta=0
    file_age=0
    if [ "$after_exists" -eq 0 ]; then
        return 0
    fi

    file_delta=$((after_size - before_size))
    file_age=$((now - after_mtime))
    [ "$file_age" -ge 0 ] || file_age=0
    if [ "$before_exists" -eq 0 ] || [ "$after_size" -ne "$before_size" ] || \
        [ "$after_mtime" -ne "$before_mtime" ]; then
        file_state="ACTIVE"
    elif [ "$file_age" -ge "$stale_after" ]; then
        file_state="STALE"
    else
        file_state="QUIET"
    fi
}

print_file_activity() {
    local label="$1"
    local path_file="$2"
    local state="$3"
    local size="$4"
    local delta="$5"
    local age="$6"

    if ! is_monitorable_path "$path_file"; then
        printf '      %-7s %s\n' "$label:" "not observable (${path_file:-not configured})"
        return 0
    fi

    printf '      %-7s %s\n' "$label:" "$path_file"
    case "$state" in
        ACTIVE)
            printf '               %s✅ ACTIVE%s  size=%s B  delta=%+d B  age=' \
                "$ok_color" "$reset" "$size" "$delta"
            format_seconds "$age"
            printf '\n'
            ;;
        STALE)
            printf '               %s⚠️  STALE%s   size=%s B  unchanged; age=' \
                "$warn_color" "$reset" "$size"
            format_seconds "$age"
            printf '\n'
            ;;
        QUIET)
            printf '               %sQUIET%s      size=%s B  unchanged; age=' \
                "$desc_color" "$reset" "$size"
            format_seconds "$age"
            printf '\n'
            ;;
        MISSING)
            printf '               %s⚠️  MISSING%s no file exists yet\n' "$warn_color" "$reset"
            ;;
    esac
}

slurm_monitor() {
    local squeue_out squeue_after job_info job_id job_name elapsed start_time partition
    local nodes cpus node_list path_workdir path_stdout_raw path_stderr_raw
    local path_stdout path_stderr snapshot now n_jobs index elapsed_seconds
    local out_state err_state job_state latest_mtime latest_age observable
    local second_snapshot_ok=1
    local n_running=0
    local n_active_job=0
    local n_quiet_job=0
    local n_suspicious_job=0
    local n_unknown_job=0
    local n_ended_job=0
    local ljob_id=() ljob_name=() lelapsed=() lstart_time=() lpartition=() lnodes=() lcpus=()
    local lnode_list=() lpath_workdir=() lpath_stdout=() lpath_stderr=()
    local lout_before_exists=() lout_before_size=() lout_before_mtime=()
    local lerr_before_exists=() lerr_before_size=() lerr_before_mtime=()
    local out_after_exists out_after_size out_after_mtime
    local err_after_exists err_after_size err_after_mtime
    local out_age out_delta err_age err_delta
    declare -A dict_still_running=()

    section "Slurm / --monitor"
    if ! squeue_out="$(squeue -h -u "$USER" -t RUNNING \
        -o '%i|%j|%M|%S|%P|%D|%C|%R' 2>&1)"; then
        fail "squeue failed while reading RUNNING jobs: $squeue_out"
    fi

    if [ -z "$squeue_out" ]; then
        printf '  %s✅ No RUNNING Slurm jobs for user %s.%s\n' "$ok_color" "$USER" "$reset"
        return 0
    fi

    while IFS='|' read -r job_id job_name elapsed start_time partition nodes cpus node_list; do
        [ -n "$job_id" ] || continue
        job_info="$(scontrol show job "$job_id" -o 2>/dev/null)"
        path_workdir="$(get_job_field "$job_info" "WorkDir" || true)"
        path_stdout_raw="$(get_job_field "$job_info" "StdOut" || true)"
        path_stderr_raw="$(get_job_field "$job_info" "StdErr" || true)"
        path_stdout="$(resolve_job_path "$path_stdout_raw" "$path_workdir")"
        path_stderr="$(resolve_job_path "$path_stderr_raw" "$path_workdir")"

        ljob_id+=("$job_id")
        ljob_name+=("$job_name")
        lelapsed+=("$elapsed")
        lstart_time+=("$start_time")
        lpartition+=("$partition")
        lnodes+=("$nodes")
        lcpus+=("$cpus")
        lnode_list+=("$node_list")
        lpath_workdir+=("$path_workdir")
        lpath_stdout+=("$path_stdout")
        lpath_stderr+=("$path_stderr")

        snapshot="$(snapshot_file "$path_stdout")"
        read -r out_after_exists out_after_size out_after_mtime <<< "$snapshot"
        lout_before_exists+=("$out_after_exists")
        lout_before_size+=("$out_after_size")
        lout_before_mtime+=("$out_after_mtime")

        snapshot="$(snapshot_file "$path_stderr")"
        read -r err_after_exists err_after_size err_after_mtime <<< "$snapshot"
        lerr_before_exists+=("$err_after_exists")
        lerr_before_size+=("$err_after_size")
        lerr_before_mtime+=("$err_after_mtime")
    done <<< "$squeue_out"

    n_jobs="${#ljob_id[@]}"
    printf '  Snapshot: %s\n' "$(date '+%F %T %Z')"
    printf '  User: %s   RUNNING at first snapshot: %d\n' "$USER" "$n_jobs"
    printf '  Sampling StdOut/StdErr for %ss; suspicious threshold: ' "$monitor_interval"
    format_seconds "$stale_after"
    printf '.\n'
    sleep "$monitor_interval"

    if ! squeue_after="$(squeue -h -u "$USER" -t RUNNING -o '%i' 2>&1)"; then
        second_snapshot_ok=0
        warn "second squeue snapshot failed; job end detection is unavailable: $squeue_after"
    else
        while IFS= read -r job_id; do
            [ -n "$job_id" ] && dict_still_running["$job_id"]=1
        done <<< "$squeue_after"
    fi

    now="$(date +%s)"
    for ((index = 0; index < n_jobs; index++)); do
        job_id="${ljob_id[$index]}"
        snapshot="$(snapshot_file "${lpath_stdout[$index]}")"
        read -r out_after_exists out_after_size out_after_mtime <<< "$snapshot"
        classify_file_activity \
            "${lout_before_exists[$index]}" "${lout_before_size[$index]}" \
            "${lout_before_mtime[$index]}" "$out_after_exists" "$out_after_size" \
            "$out_after_mtime" "$now" "$stale_after"
        out_state="$file_state"
        out_delta="$file_delta"
        out_age="$file_age"

        if [ "${lpath_stderr[$index]}" = "${lpath_stdout[$index]}" ]; then
            err_state="SAME"
            err_after_exists=0
            err_after_size=0
            err_after_mtime=0
            err_delta=0
            err_age=0
        else
            snapshot="$(snapshot_file "${lpath_stderr[$index]}")"
            read -r err_after_exists err_after_size err_after_mtime <<< "$snapshot"
            classify_file_activity \
                "${lerr_before_exists[$index]}" "${lerr_before_size[$index]}" \
                "${lerr_before_mtime[$index]}" "$err_after_exists" "$err_after_size" \
                "$err_after_mtime" "$now" "$stale_after"
            err_state="$file_state"
            err_delta="$file_delta"
            err_age="$file_age"
        fi

        observable=0
        is_monitorable_path "${lpath_stdout[$index]}" && observable=1
        if [ "$err_state" != "SAME" ] && is_monitorable_path "${lpath_stderr[$index]}"; then
            observable=1
        fi
        latest_mtime=0
        [ "$out_after_exists" -eq 1 ] && latest_mtime="$out_after_mtime"
        if [ "$err_after_exists" -eq 1 ] && [ "$err_after_mtime" -gt "$latest_mtime" ]; then
            latest_mtime="$err_after_mtime"
        fi
        latest_age=$((now - latest_mtime))
        elapsed_seconds="$(duration_to_seconds "${lelapsed[$index]}")"

        if [ "$second_snapshot_ok" -eq 1 ] && [ -z "${dict_still_running[$job_id]:-}" ]; then
            job_state="ENDED"
            n_ended_job=$((n_ended_job + 1))
        elif [ "$out_state" = "ACTIVE" ] || [ "$err_state" = "ACTIVE" ]; then
            job_state="ACTIVE"
            n_active_job=$((n_active_job + 1))
            n_running=$((n_running + 1))
        elif [ "$observable" -eq 0 ]; then
            job_state="UNKNOWN"
            n_unknown_job=$((n_unknown_job + 1))
            n_running=$((n_running + 1))
        elif { [ "$latest_mtime" -gt 0 ] && [ "$latest_age" -ge "$stale_after" ]; } || \
            { [ "$latest_mtime" -eq 0 ] && [ "$elapsed_seconds" -ge "$stale_after" ]; }; then
            job_state="SUSPICIOUS"
            n_suspicious_job=$((n_suspicious_job + 1))
            n_running=$((n_running + 1))
        else
            job_state="QUIET"
            n_quiet_job=$((n_quiet_job + 1))
            n_running=$((n_running + 1))
        fi

        printf '\n  %sJobID=%s%s  Name=%s  Status=' "$bold" "$job_id" "$reset" "${ljob_name[$index]}"
        case "$job_state" in
            ACTIVE) printf '%s✅ ACTIVE%s\n' "$ok_color" "$reset" ;;
            SUSPICIOUS) printf '%s⚠️  SUSPICIOUS%s\n' "$warn_color" "$reset" ;;
            ENDED) printf '%sENDED during sample%s\n' "$desc_color" "$reset" ;;
            UNKNOWN) printf '%sUNKNOWN%s\n' "$warn_color" "$reset" ;;
            *) printf '%sQUIET%s\n' "$desc_color" "$reset" ;;
        esac
        printf '      Elapsed: %-12s Partition: %-12s Nodes: %s  CPUs: %s  NodeList: %s\n' \
            "${lelapsed[$index]}" "${lpartition[$index]}" "${lnodes[$index]}" \
            "${lcpus[$index]}" "${lnode_list[$index]}"
        printf '      Start:   %s\n' "${lstart_time[$index]}"
        printf '      WorkDir: %s\n' "${lpath_workdir[$index]:-unknown}"
        print_file_activity "StdOut" "${lpath_stdout[$index]}" "$out_state" \
            "$out_after_size" "$out_delta" "$out_age"
        if [ "$err_state" = "SAME" ]; then
            printf '      %-7s same as StdOut\n' "StdErr:"
        else
            print_file_activity "StdErr" "${lpath_stderr[$index]}" "$err_state" \
                "$err_after_size" "$err_delta" "$err_age"
        fi
    done

    printf '\n%s================ 📊 Slurm monitor summary%s\n' "$bold" "$reset"
    printf '  Running now: %d   %s✅ active: %d%s   quiet: %d   %s⚠️  suspicious: %d%s' \
        "$n_running" "$ok_color" "$n_active_job" "$reset" "$n_quiet_job" \
        "$warn_color" "$n_suspicious_job" "$reset"
    printf '   unknown: %d   ended: %d\n' "$n_unknown_job" "$n_ended_job"
    printf '  %sSUSPICIOUS means no observed Slurm output within the threshold; buffering or a long compute step can be normal.%s\n' \
        "$desc_color" "$reset"
}
# to here

################ cheat sheet ################

print_commands() {
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
    command_item "$script_name --monitor" \
        "Summarize RUNNING jobs and sample StdOut/StdErr activity."

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
    command_item "$script_name --openai" \
        "Probe the OpenAI endpoints (Codex / ChatGPT) from this node."
    command_item "$script_name --claude" \
        "Probe the Anthropic endpoints (Claude Code / claude.ai) from this node."
    command_item "$script_name --net --ip --proxy" \
        "Probe every endpoint group, show the public IP and the proxy variables."
    command_item "$script_name --help" \
        "List all options of this script."
}
# to here

################ usage ################

usage() {
    printf '%sUsage: %s [options]%s\n\n' "$bold" "$script_name" "$reset"
    printf 'Without options, print the useful-command cheat sheet.\n'
    printf 'With options, probe the network, monitor Slurm jobs, or print selected information.\n\n'
    printf 'Options:\n'
    printf '  -h, --help         Show this help and exit.\n'
    printf '      --openai       Probe OpenAI endpoints (api.openai.com, chatgpt.com).\n'
    printf '      --claude       Probe Anthropic endpoints (api.anthropic.com, claude.ai).\n'
    printf '      --github       Probe GitHub web / API / raw download.\n'
    printf '      --pypi         Probe PyPI, its CDN and the Tsinghua mirror.\n'
    printf '      --conda        Probe Anaconda / conda-forge and the Tsinghua mirror.\n'
    printf '      --base         Probe baseline sites (baidu / google / 1.1.1.1).\n'
    printf '      --net          Probe all groups above.\n'
    printf '      --ip           Show public IP information (ipinfo.io, cip.cc).\n'
    printf '      --proxy        Show the proxy environment variables.\n'
    printf '      --monitor      Show RUNNING jobs, paths, resources and output activity.\n'
    printf '      --monitor-interval SEC\n'
    printf '                      Seconds between output-file samples (default: %s).\n' "$monitor_interval_default"
    printf '      --stale-after SEC\n'
    printf '                      Warn after this many seconds without output (default: %s).\n' "$stale_after_default"
    printf '      --commands     Print the cheat sheet (implicit when no option is given).\n'
    printf '      --all          --net --ip --proxy --commands.\n'
    printf '      --timeout SEC  Per-request timeout in seconds (default: %s).\n' "$timeout_default"
}
# to here

################ check: parse and validate the arguments ################

timeout="$timeout_default"
monitor_interval="$monitor_interval_default"
stale_after="$stale_after_default"
lgroup_selected=()
show_ip=0
show_proxy=0
show_commands=0
show_monitor=0

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --base|--openai|--claude|--github|--pypi|--conda)
            lgroup_selected+=("${1#--}")
            ;;
        --net)
            lgroup_selected+=("${lgroup_order[@]}")
            ;;
        --ip)
            show_ip=1
            ;;
        --proxy)
            show_proxy=1
            ;;
        --monitor)
            show_monitor=1
            ;;
        --monitor-interval)
            [ $# -ge 2 ] || fail "--monitor-interval needs a value in seconds"
            monitor_interval="$(check_positive_int "$2" "--monitor-interval")" || exit 1
            shift
            ;;
        --monitor-interval=*)
            monitor_interval="$(check_positive_int "${1#--monitor-interval=}" "--monitor-interval")" || exit 1
            ;;
        --stale-after)
            [ $# -ge 2 ] || fail "--stale-after needs a value in seconds"
            stale_after="$(check_positive_int "$2" "--stale-after")" || exit 1
            shift
            ;;
        --stale-after=*)
            stale_after="$(check_positive_int "${1#--stale-after=}" "--stale-after")" || exit 1
            ;;
        --commands)
            show_commands=1
            ;;
        --all)
            lgroup_selected+=("${lgroup_order[@]}")
            show_ip=1
            show_proxy=1
            show_commands=1
            ;;
        --timeout)
            [ $# -ge 2 ] || fail "--timeout needs a value in seconds"
            timeout="$(check_positive_int "$2" "--timeout")" || exit 1
            shift
            ;;
        --timeout=*)
            timeout="$(check_positive_int "${1#--timeout=}" "--timeout")" || exit 1
            ;;
        *)
            printf '%s❌ ERROR: unknown option: %s%s\n\n' "$fail_color" "$1" "$reset" >&2
            usage >&2
            exit 1
            ;;
    esac
    shift
done

# No option at all keeps the historical behaviour: just print the cheat sheet.
if [ "${#lgroup_selected[@]}" -eq 0 ] && [ "$show_ip" -eq 0 ] && \
    [ "$show_proxy" -eq 0 ] && [ "$show_monitor" -eq 0 ]; then
    show_commands=1
fi

if [ "${#lgroup_selected[@]}" -gt 0 ] || [ "$show_ip" -eq 1 ]; then
    command -v curl >/dev/null 2>&1 || fail "curl is not available, cannot probe the network"
fi
if [ "$show_monitor" -eq 1 ]; then
    command -v squeue >/dev/null 2>&1 || fail "squeue is not available, cannot monitor Slurm jobs"
    command -v scontrol >/dev/null 2>&1 || fail "scontrol is not available, cannot inspect Slurm job paths"
    command -v stat >/dev/null 2>&1 || fail "stat is not available, cannot inspect Slurm output files"
fi
# to here

################ main ################

[ "$show_commands" -eq 1 ] && print_commands
[ "$show_monitor" -eq 1 ] && slurm_monitor

if [ "${#lgroup_selected[@]}" -gt 0 ]; then
    # Deduplicate while keeping the canonical order (--net --claude must not probe twice).
    for group in "${lgroup_order[@]}"; do
        for selected in "${lgroup_selected[@]}"; do
            if [ "$group" = "$selected" ]; then
                net_check_group "$group"
                break
            fi
        done
    done
fi

[ "$show_proxy" -eq 1 ] && net_show_proxy
[ "$show_ip" -eq 1 ] && net_show_ip

if [ "${#lgroup_selected[@]}" -gt 0 ]; then
    print_summary || exit 1
fi

exit 0
# to here
