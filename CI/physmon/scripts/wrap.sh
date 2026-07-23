#!/usr/bin/env bash
#
# wrap.sh --witness <path> --step <name> [--log <path>] [--touches "<list>"] -- <cmd...>
#
# Runs <cmd>; writes a small JSON witness {"step", "rc", "log"} to <witness>;
# touches any path in <list> that does not already exist (so snakemake's
# declared outputs are always present and the DAG can continue even when an
# upstream step failed); always exits 0 itself.
#
# Used by every rule in the physmon Snakefile so the workflow runs to
# completion regardless of per-step failures; failure surfacing happens in the
# summary step which reads the witness JSON files.

set -u

# Restore the full ACTS runtime environment that run_snakemake.sh moved aside so
# the snakemake driver could start (see scripts/run_snakemake.sh). The wrapped
# generator/histcmp/plot command needs ROOT/dd4hep/acts on LD_LIBRARY_PATH and
# their Python packages on PYTHONPATH. Absent (a bare `snakemake` invocation) →
# leave the inherited environment untouched.
if [ -n "${_PHYSMON_JOB_LD_LIBRARY_PATH+x}" ]; then
    export LD_LIBRARY_PATH="$_PHYSMON_JOB_LD_LIBRARY_PATH"
fi
if [ -n "${_PHYSMON_JOB_PYTHONPATH+x}" ]; then
    export PYTHONPATH="$_PHYSMON_JOB_PYTHONPATH"
fi

witness=""
step=""
log_path=""
touches=""

while [ $# -gt 0 ]; do
    case "$1" in
        --witness) witness="$2"; shift 2;;
        --step)    step="$2";    shift 2;;
        --log)     log_path="$2"; shift 2;;
        --touches) touches="$2"; shift 2;;
        --)        shift; break;;
        *) echo "wrap.sh: unknown arg: $1" >&2; exit 2;;
    esac
done

if [ -z "$witness" ] || [ -z "$step" ] || [ $# -eq 0 ]; then
    echo "wrap.sh: usage: $0 --witness <p> --step <n> [--log <p>] [--touches '<list>'] -- <cmd...>" >&2
    exit 2
fi

mkdir -p "$(dirname "$witness")"

# Run the command without aborting on non-zero. If a log path was given, redirect
# the command's stdout/stderr there. wrap.sh's own messages go to its inherited
# stderr (snakemake's main log).
set +e
if [ -n "$log_path" ]; then
    mkdir -p "$(dirname "$log_path")"
    : > "$log_path"
    "$@" >> "$log_path" 2>&1
else
    "$@"
fi
rc=$?

# JSON-escape the log path (only need to handle backslash and double quote;
# we control the character set).
if [ -n "$log_path" ]; then
    escaped=${log_path//\\/\\\\}
    escaped=${escaped//\"/\\\"}
    log_json="\"$escaped\""
else
    log_json="null"
fi

printf '{"step":"%s","rc":%d,"log":%s}\n' "$step" "$rc" "$log_json" > "$witness"

# Touch any declared outputs that didn't get produced. The witness itself was
# just written above; skip it. Empty 'touches' is fine — the for loop is a no-op.
for f in $touches; do
    [ "$f" = "$witness" ] && continue
    if [ ! -e "$f" ]; then
        mkdir -p "$(dirname "$f")"
        touch "$f"
    fi
done

exit 0
