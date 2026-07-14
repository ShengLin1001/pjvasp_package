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
#   pei_slurm_univ_useful_command.sh --all --timeout 3

script_name="$(basename "$0")"
timeout_default=8

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
    printf 'With options, probe the network instead (exit 1 if any endpoint is unreachable).\n\n'
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
    printf '      --commands     Print the cheat sheet (implicit when no option is given).\n'
    printf '      --all          --net --ip --proxy --commands.\n'
    printf '      --timeout SEC  Per-request timeout in seconds (default: %s).\n' "$timeout_default"
}
# to here

################ check: parse and validate the arguments ################

timeout="$timeout_default"
lgroup_selected=()
show_ip=0
show_proxy=0
show_commands=0

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
if [ "${#lgroup_selected[@]}" -eq 0 ] && [ "$show_ip" -eq 0 ] && [ "$show_proxy" -eq 0 ]; then
    show_commands=1
fi

if [ "${#lgroup_selected[@]}" -gt 0 ] || [ "$show_ip" -eq 1 ]; then
    command -v curl >/dev/null 2>&1 || fail "curl is not available, cannot probe the network"
fi
# to here

################ main ################

[ "$show_commands" -eq 1 ] && print_commands

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
