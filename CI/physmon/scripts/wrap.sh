#!/usr/bin/env bash
#
# wrap.sh --witness <path> --step <name> [--log <path>] [--touches "<paths>"] -- <cmd...>
#
# Runs <cmd>, with its output in <log> if given, and writes the witness
# {"step", "exit_code", "log"}. Creates every path in <paths> that <cmd> did
# not, and always exits 0. The physmon rules run through this, so the workflow
# carries on after a failed step, and the summary reports the failure from the
# witnesses.

set -u

witness="" step="" log="" touches=""
while [ $# -gt 0 ]; do
    case "$1" in
        --witness) witness="$2"; shift 2;;
        --step)    step="$2";    shift 2;;
        --log)     log="$2";     shift 2;;
        --touches) touches="$2"; shift 2;;
        --)        shift; break;;
        *) echo "wrap.sh: unknown argument: $1" >&2; exit 2;;
    esac
done
if [ -z "$witness" ] || [ -z "$step" ] || [ $# -eq 0 ]; then
    echo "wrap.sh: usage: $0 --witness <path> --step <name> [--log <path>] [--touches '<paths>'] -- <cmd...>" >&2
    exit 2
fi

if [ -n "$log" ]; then
    mkdir -p "$(dirname "$log")"
    "$@" > "$log" 2>&1
    exit_code=$?
    log_json="\"$log\""
else
    "$@"
    exit_code=$?
    log_json="null"
fi

mkdir -p "$(dirname "$witness")"
printf '{"step":"%s","exit_code":%d,"log":%s}\n' "$step" "$exit_code" "$log_json" > "$witness"

for path in $touches; do
    if [ ! -e "$path" ]; then
        mkdir -p "$(dirname "$path")"
        touch "$path"
    fi
done

exit 0
