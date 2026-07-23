#!/usr/bin/env bash
# Checks that every checked-in requirements.txt lockfile can actually be
# installed on every Python version Acts supports.
#
# The lockfiles are compiled with `uv pip compile --universal`, which keeps
# environment markers so one file covers the whole supported range. This check
# guards that property: without it, a lockfile compiled for a single version
# silently breaks every other interpreter, which is only noticed by the nightly
# PyPI build (PR CI builds a single Python version).
#
# Usage:
#   CI/check_requirements.sh              # check all supported versions
#   CI/check_requirements.sh 3.11 3.14    # check specific versions
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"

# Keep in sync with the folder_list in
# .github/workflows/update-pip-requirements.yml
lockfiles=(
  CI/clang_tidy/requirements.txt
  CI/fpe_masks/requirements.txt
  codegen/requirements.txt
  docs/old/requirements.txt
  Python/Examples/tests/requirements.txt
  Examples/Scripts/requirements.txt
)

versions=("$@")
if [[ ${#versions[@]} -eq 0 ]]; then
  # Derive the version list from `requires-python` in pyproject.toml so the
  # check stays in sync when that specifier changes.
  mapfile -t versions < <("${SCRIPT_DIR}/python_versions.py" --list)
fi

echo "Supported Python versions: ${versions[*]}"

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT

ec=0
for version in "${versions[@]}"; do
  venv="${workdir}/venv-${version}"
  echo "::group::Python ${version}"
  uv venv --seed --python "${version}" "${venv}"

  for lockfile in "${lockfiles[@]}"; do
    # A dry run resolves and checks wheel availability without downloading
    # anything, which is all we need to catch version-incompatible pins.
    if "${venv}/bin/python" -m pip install --dry-run --quiet \
      -r "${REPO_ROOT}/${lockfile}"; then
      echo "PASS ${lockfile} on Python ${version}"
    else
      echo "FAIL ${lockfile} on Python ${version}"
      ec=1
    fi
  done
  echo "::endgroup::"
done

if [[ ${ec} -ne 0 ]]; then
  echo
  echo "Some lockfiles are not installable on all supported Python versions."
  echo "Regenerate them with:"
  echo "  uv pip compile --universal --python-version \$(CI/python_versions.py --floor) ..."
  echo "see .github/workflows/update-pip-requirements.yml"
fi

exit ${ec}
