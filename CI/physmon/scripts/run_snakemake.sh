#!/usr/bin/env bash
#
# run_snakemake.sh [snakemake args...]
#
# Launch the snakemake *driver* with an environment that does NOT contain the
# ACTS runtime libraries, then let the individual jobs restore that environment.
#
# Why: physmon runs after `source build/this_acts_withdeps.sh`, which prepends
# ROOT/dd4hep/podio to LD_LIBRARY_PATH (and their Python packages to PYTHONPATH).
# Those are built against the ACTS Python and are ABI-incompatible with the
# interpreter uv installs snakemake into, so the snakemake process segfaults
# (exit 139) at startup, before any rule runs. Neutralising PYTHONPATH alone
# (python -E) is not enough — the crash is at the dynamic-linker level via
# LD_LIBRARY_PATH.
#
# The CI step snapshots the pristine env (before sourcing the ACTS env) into
# _PHYSMON_DRIVER_*; we run the driver with that. The full ACTS env is passed
# through to the jobs snakemake launches via _PHYSMON_JOB_* — scripts/wrap.sh
# (generators/histcmp/plot) and the summary rule re-export it. If the caller
# provided no snapshot (e.g. a local run that doesn't crash), we fall back to
# the current environment, with `-E` still shielding the interpreter from the
# ACTS PYTHONPATH.
set -euo pipefail

# Full ACTS runtime env → handed to the jobs (see wrap.sh / the summary rule).
# The workflow scripts (run from workflows/) import helper modules that live
# outside the ACTS Python package: physmon_common from CI/physmon, and the
# truth_tracking_* example drivers from Examples/Scripts/Python. Append both so
# those imports resolve — appended last so they can never shadow the ACTS/ROOT
# packages earlier on the path. Relative paths resolve since the jobs run from
# the repo root (snakemake's working directory).
export _PHYSMON_JOB_LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
export _PHYSMON_JOB_PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}CI/physmon:Examples/Scripts/Python"

snakemake_bin=$(command -v snakemake)
# Interpreter from snakemake's shebang, so this works wherever it was installed
# (e.g. `uv tool install`) without hard-coding a path.
snakemake_py=$(sed -n '1s/^#!//p' "$snakemake_bin")

exec env \
  LD_LIBRARY_PATH="${_PHYSMON_DRIVER_LD_LIBRARY_PATH-${LD_LIBRARY_PATH:-}}" \
  PYTHONPATH="${_PHYSMON_DRIVER_PYTHONPATH-${PYTHONPATH:-}}" \
  "$snakemake_py" -E "$snakemake_bin" "$@"
