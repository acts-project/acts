# Environment of a physmon rule's shell, sourced by the Snakefile for every
# rule: the ACTS build in <repo>/build, plus the two directories the workflow
# scripts import from. Activating it per rule keeps the ACTS PYTHONPATH out of
# the snakemake process, whose interpreter would crash on its ROOT, dd4hep and
# podio modules.

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REPO_ROOT=$( cd -- "${SCRIPT_DIR}/../.." &> /dev/null && pwd )

if [ ! -f "${REPO_ROOT}/build/this_acts_withdeps.sh" ]; then
    echo "physmon: no ACTS build at ${REPO_ROOT}/build (this_acts_withdeps.sh not found)" >&2
    return 1
fi

# this_acts_withdeps.sh reads unset variables, which the rule shells' `set -u`
# forbids: switch it off for the source and restore it as it was.
_physmon_nounset=$(set +o | grep nounset)
set +u
source "${REPO_ROOT}/build/this_acts_withdeps.sh"
eval "${_physmon_nounset}"
unset _physmon_nounset

# Appended, so they cannot shadow ACTS or ROOT.
export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}${REPO_ROOT}/CI/physmon:${REPO_ROOT}/Examples/Scripts/Python"
