# Environment for the physmon rules, stacked on top of a plain shell:
#
#   1. the ACTS build environment (build/this_acts_withdeps.sh), which is what
#      puts acts, ROOT, dd4hep and podio on PYTHONPATH,
#   2. the two import roots the workflow scripts need that are not part of the
#      ACTS Python package: physmon_common from CI/physmon, and the
#      truth_tracking_* drivers from Examples/Scripts/Python.
#
# The Snakefile sources this from every rule's shell rather than expecting the
# caller to have activated anything first. That keeps the ACTS PYTHONPATH out of
# the snakemake process itself: the ROOT/dd4hep/podio extension modules on it are
# ABI-tied to the ACTS interpreter and segfault any other interpreter that
# imports them, which used to mean snakemake had to be installed against the
# ACTS Python.
#
# The build directory is assumed to be <repo>/build.

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPO_ROOT=$( cd -- "${SCRIPT_DIR}/../.." &> /dev/null && pwd )
BUILD_DIR="${REPO_ROOT}/build"

if [ ! -f "${BUILD_DIR}/this_acts_withdeps.sh" ]; then
    echo "physmon: no ACTS build environment at ${BUILD_DIR}/this_acts_withdeps.sh" >&2
    echo "physmon: the workflow expects the build directory at ${BUILD_DIR}" >&2
    return 1
fi

# this_acts_withdeps.sh expands ${LD_LIBRARY_PATH} and friends unguarded, so it
# cannot run under `set -u` — which the rule shells do set. Drop it for the
# source and put it back exactly as it was.
case $- in
    *u*) _physmon_nounset=1 ;;
    *)   _physmon_nounset=0 ;;
esac
set +u
source "${BUILD_DIR}/this_acts_withdeps.sh"
if [ "${_physmon_nounset}" = "1" ]; then
    set -u
fi
unset _physmon_nounset

# Appended, so they cannot shadow ACTS/ROOT.
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}${REPO_ROOT}/CI/physmon:${REPO_ROOT}/Examples/Scripts/Python"
