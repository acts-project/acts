#!/usr/bin/env bash
# Builds ACTS Python wheels via cibuildwheel with the project's standard
# CIBW_* configuration. All arguments are forwarded to cibuildwheel
# unchanged (e.g. a package-dir positional, --only, --output-dir, ...).
# Anything else is configured via environment variables: CIBW_BUILD selects
# which target(s) to build and must be set by the caller; GITHUB_TOKEN is
# forwarded into the build environment if set.
set -euo pipefail

if [[ -z "${CIBW_BUILD:-}" ]]; then
  echo "CIBW_BUILD must be set" >&2
  exit 1
fi

CI="${CI:-true}"
export CI
# Resolved on the host so both platforms agree on where the cache lives.
CCACHE_DIR="${CCACHE_DIR:-$PWD/ccache}"
export CCACHE_DIR
export CIBW_MANYLINUX_X86_64_IMAGE="manylinux_2_34" # based on almalinux9
export CIBW_SKIP="*-musllinux* *-manylinux_i686"
SETUP_CMD="bash {package}/CI/dependencies/setup.sh -t v23.3.1 -d deps -e env.sh"
export CIBW_BEFORE_ALL_LINUX="dnf install -y bc ccache && ${SETUP_CMD}"
export CIBW_BEFORE_ALL_MACOS="brew install ninja ccache && ${SETUP_CMD}"
export CIBW_ENVIRONMENT="GITHUB_TOKEN=${GITHUB_TOKEN:-}"
export CIBW_ENVIRONMENT_PASS="CI"
export CIBW_BEFORE_BUILD="ccache -z"
export CIBW_ENVIRONMENT_LINUX="CMAKE_PREFIX_PATH=\$PWD/deps/venv:\$PWD/deps/view CCACHE_DIR=/host${CCACHE_DIR} LD_LIBRARY_PATH=\$PWD/deps/view/lib64:\$PWD/deps/view/lib:\$PWD/deps/venv/lib64:\$PWD/deps/venv/lib"
export CIBW_ENVIRONMENT_MACOS="CMAKE_PREFIX_PATH=\$PWD/deps/venv:\$PWD/deps/view CCACHE_DIR=${CCACHE_DIR} MACOSX_DEPLOYMENT_TARGET=26.0"
export CIBW_BEFORE_TEST="ccache -s && uv pip install -r {package}/Python/Examples/tests/requirements.txt"
export CIBW_TEST_COMMAND="pytest {package}/Python/Examples/tests -m pypi -v"
# patchelf 0.17.2 (pinned in the manylinux image) corrupts auditwheel-vendored
# libs (e.g. libzstd) it repairs, causing a segfault at import time. Force a
# newer patchelf until manylinux ships a stable release with the fix.
export CIBW_REPAIR_WHEEL_COMMAND_LINUX="pipx install --force --pip-args='--pre' patchelf==0.19.0.0rc1 && auditwheel repair -w {dest_dir} {wheel}"

uv tool run cibuildwheel==3.4.1 "$@"
