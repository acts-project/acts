#!/bin/bash
set -u
set -e
set -o pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

export SPACK_COLOR=always

function start_section() {
    local section_name="$1"
    if [ -n "${GITHUB_ACTIONS:-}" ]; then
        echo "::group::${section_name}"
    else
        echo "${section_name}"
    fi
}

function end_section() {
    if [ -n "${GITHUB_ACTIONS:-}" ]; then
        echo "::endgroup::"
    fi
}

# Signatures in a command's output that indicate a *transient* failure (network
# / registry hiccup) that a retry can plausibly fix. Anything that does not
# match one of these is treated as a genuine failure and is NOT retried, so we
# don't burn CI time looping on an error a retry cannot resolve (compile error,
# disk full, bad spec, ...).
#
# NOTE on "no binary available": spack does *not* surface a throttled/failed
# GHCR fetch as a network error. Under `--use-buildcache only` it swallows the
# fetch failure and reports the spec as simply having no binary, which is by far
# the most common transient failure we see. So those signatures are treated as
# transient and retried. The cost of being wrong is bounded (a genuinely missing
# binary just burns the retry budget before failing), whereas not retrying turns
# every GHCR hiccup into a red build.
TRANSIENT_ERROR_PATTERNS=(
  # spack buildcache misses caused by upstream fetch/rate-limit failures
  "No binary for"
  "no binary available"
  "NoBinaryFoundError"
  "cannot install.*cache-only"
  "Failed to install.*due to.*FetchError"
  "Failed to fetch"
  "fetch failed"
  "FetchError"
  "curl.*(error|failed|timed out|Could not resolve)"
  "Connection reset"
  "Connection refused"
  "Connection timed out"
  "Could not resolve host"
  "Temporary failure in name resolution"
  "Network is unreachable"
  "Read timed out"
  "read timeout"
  "timed out"
  "[Tt]imeout"
  "TLS handshake"
  "EOF occurred"
  "Server disconnected"
  "Remote end closed connection"
  "toomanyrequests"
  "Too Many Requests"
  "rate limit"
  "HTTP Error 429"
  "HTTP Error 5[0-9][0-9]"
  "50[0-9] (Internal Server Error|Bad Gateway|Service Unavailable|Gateway Time-out)"
  "unable to access.*github.com"
  "RPC failed"
  "early EOF"
)

# Run "$@", streaming its combined output while also capturing it. On failure,
# retry with exponential backoff *only* when the captured output matches a known
# transient-error signature; otherwise return the command's exit status
# immediately. Output is streamed live (via tee) so long-running commands don't
# trip CI "no output" watchdogs.
#
# Backoff is 20s, 40s, 80s, 160s, 320s (+/- jitter) by default: ~10 min of total
# wait across 6 attempts. GHCR rate-limit windows outlast a short backoff, so a
# tight schedule just burns all attempts inside the same bad window and still
# fails. Jitter keeps the many parallel CI jobs from retrying in lockstep and
# re-triggering the limit together.
function retry_transient() {
  local max_attempts=${DEP_RETRY_MAX_ATTEMPTS:-6}
  local delay=${DEP_RETRY_BASE_DELAY:-20}
  local attempt=1
  local log status pat matched jitter sleep_for
  log=$(mktemp)
  while true; do
    echo "attempt ${attempt}/${max_attempts}: $*"
    status=0
    # pipefail (set above) makes PIPESTATUS[0] the command's own exit status.
    "$@" 2>&1 | tee "${log}" || status=${PIPESTATUS[0]}
    if [ "${status}" -eq 0 ]; then
      rm -f "${log}"
      return 0
    fi

    matched=""
    for pat in "${TRANSIENT_ERROR_PATTERNS[@]}"; do
      if grep -qiE -- "${pat}" "${log}"; then
        matched="${pat}"
        break
      fi
    done
    rm -f "${log}"

    if [ -z "${matched}" ]; then
      echo "Command failed (exit ${status}) with no transient-error signature; not retrying"
      return "${status}"
    fi
    if [ "${attempt}" -ge "${max_attempts}" ]; then
      echo "Command still failing after ${max_attempts} attempts (last transient signature: '${matched}')"
      return "${status}"
    fi
    # +/-25% jitter so parallel jobs don't retry in lockstep.
    jitter=$(( (RANDOM % (delay / 2 + 1)) - delay / 4 ))
    sleep_for=$(( delay + jitter ))
    [ "${sleep_for}" -lt 1 ] && sleep_for=1
    echo "Transient failure detected (matched '${matched}'); retrying in ${sleep_for}s"
    sleep "${sleep_for}"
    attempt=$((attempt + 1))
    delay=$((delay * 2))
  done
}

# Parse command line arguments
while getopts "c:t:d:e:s:F:fh" opt; do
  case ${opt} in
    c )
      compiler=$OPTARG
      ;;
    F )
      flavor=$OPTARG
      ;;
    t )
      tag=$OPTARG
      ;;
    d )
      destination=$OPTARG
      ;;
    e )
      env_file=$OPTARG
      ;;
    s )
      cxx_std=$OPTARG
      ;;
    f )
      full_install=true
      ;;
    h )
      echo "Usage: $0 [-c compiler] [-t tag] [-d destination] -e env_file [-h]"
      echo "Options:"
      echo "  -c <compiler>    Specify compiler (defaults to CXX env var)"
      echo "  -t <tag>         Specify dependency tag (defaults to DEPENDENCY_TAG env var)"
      echo "  -d <destination> Specify install destination (defaults based on CI environment)"
      echo "  -e <env_file>    Specify environment file to output environments to"
      echo "  -s <cxx_std>     C++ standard for lockfile selection (e.g. 20, 23). Defaults to CXXSTD env var or 20."
      echo "  -F <flavor>      Accelerator flavor (e.g. cuda13, rocm-gfx90a). Defaults to FLAVOR env var or the host stack."
      echo "  -f               Full dependency installation. Includes Geant4 datasets and Python packages."
      echo "  -h               Show this help message"
      exit 0
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument" 1>&2
      exit 1
      ;;
  esac
done

script_start=$(date +%s.%N)

# Helper to print elapsed time since previous checkpoint
checkpoint() {
    local label=$1
    local now
    now=$(date +%s.%N)
    local elapsed
    elapsed=$(echo "$now - ${last_time:-$script_start}" | bc)
    printf "[%s] %.3f s\n" "$label" "$elapsed"
    last_time=$now
}

# Set defaults if not specified
if [ -z "${compiler:-}" ]; then
  compiler="${CXX:-default}"
fi

if [ -z "${tag:-}" ]; then
  tag="${DEPENDENCY_TAG:-}"
  if [ -z "${tag:-}" ]; then
    echo "No tag specified via -t or DEPENDENCY_TAG environment variable"
    exit 1
  fi
fi

if [ -z "${destination:-}" ]; then
  if [ -n "${GITHUB_ACTIONS:-}" ]; then
    destination="${GITHUB_WORKSPACE}/dependencies"
  elif [ -n "${GITLAB_CI:-}" ]; then
    destination="${CI_PROJECT_DIR}/dependencies"
  else
    echo "No destination specified via -d and not running in CI"
    exit 1
  fi
fi

if [ -z "${env_file:-}" ]; then
  echo "No environment file specified via -e"
  exit 1
fi

if [ -z "${cxx_std:-}" ]; then
  cxx_std="${CXXSTD:-20}"
fi

# `host` is the plain CPU stack, published with no flavor token: pass nothing.
if [ -z "${flavor:-}" ]; then
  flavor="${FLAVOR:-}"
fi
if [ "${flavor}" == "host" ]; then
  flavor=""
fi

checkpoint "Create environment file $(realpath "$env_file")"
echo "" > "$env_file"
export env_file

function set_env {
  key="$1"
  value="$2"

  echo "=> ${key}=${value}"

  echo "export ${key}=${value}" >> "$env_file"
}



checkpoint "Starting setup script"

echo "Install tag: $tag"
echo "Install destination: $destination"

mkdir -p "${destination}"

if [ -n "${GITLAB_CI:-}" ]; then
    _spack_folder=${CI_PROJECT_DIR}/spack
else
    _spack_folder=${PWD}/spack
fi

start_section "Install spack if not already installed"
if ! command -v spack &> /dev/null; then
  "${SCRIPT_DIR}/setup_spack.sh" "${_spack_folder}"
  source "${_spack_folder}/share/spack/setup-env.sh"
fi
checkpoint "Spack install complete"

_spack_repo_version=${SPACK_REPO_VERSION:-develop}
_spack_repo_directory="$(realpath "$(spack location --repo builtin)/../../../")"

echo "Ensure builtin repo is synced to commit ${_spack_repo_version}"

git config --global --add safe.directory "${_spack_repo_directory}"
retry_transient spack repo update builtin --commit "${_spack_repo_version}"
checkpoint "Spack repository updated"

end_section

if [ -n "${GITLAB_CI:-}" ]; then
  # Use the project spack config for GitLab CI so we can cache it
  mkdir -p ${CI_PROJECT_DIR}/.spack
  ln -s ${CI_PROJECT_DIR}/.spack ${HOME}/.spack
fi



if [ -n "${CI:-}" ]; then
  start_section "Add buildcache mirror"
  mirror_name="acts-spack-buildcache"
  mirror_url="oci://ghcr.io/acts-project/spack-buildcache"
  if [ -n "${GITLAB_CI:-}" ]; then
  # Use CERN mirror for non-Github Actions
    mirror_url="oci://registry.cern.ch/ghcr.io/acts-project/spack-buildcache"
  fi

  # Check if this buildcache is already configured
  if ! spack mirror list | grep -q ${mirror_name}; then
    echo "Adding buildcache ${mirror_name}"
    spack mirror add ${mirror_name} ${mirror_url} --unsigned
  fi
  # Authenticate GHCR reads to avoid anonymous rate limits (which spack
  # misclassifies as "no binary available"). Idempotent on cached spack installs.
  # GITHUB_TOKEN must be in env when `spack install` later fetches from the mirror.
  if [ -n "${GITHUB_TOKEN:-}" ] && [[ "${mirror_url}" == oci://ghcr.io/* ]]; then
    echo "Setting GHCR credentials on ${mirror_name}"
    spack mirror set \
      --oci-username "${GITHUB_ACTOR:-x-access-token}" \
      --oci-password-variable GITHUB_TOKEN \
      "${mirror_name}"
  fi
  # Verify the mirror config (password-variable stores only the env var name,
  # not the secret value, so this is safe to print).
  spack mirror list
  spack config get mirrors
  end_section

  start_section "Add ACTS package repository"
  if ! spack repo list | grep -q "acts"; then
    echo "Adding ACTS package repository from ci-dependencies"
    retry_transient spack repo add https://github.com/acts-project/ci-dependencies.git --path spack_repo/acts
  fi
  echo "Updating ACTS package repository to tag ${tag}"
  retry_transient spack repo update acts --tag "${tag}"
  end_section

  start_section "Locate OpenGL"
  "${SCRIPT_DIR}/opengl.sh"
  checkpoint "OpenGL location complete"
  end_section
fi

start_section "Get spack lock file"
arch=$(spack arch --family)

env_dir="${destination}/env"
view_dir="${destination}/view"
venv_dir="${destination}/venv"
mkdir -p ${env_dir}

lock_file_path="${destination}/spack.lock"
cmd=(
    "${SCRIPT_DIR}/select_lockfile.py"
    "--tag" "${tag}"
    "--arch" "${arch}"
    "--cxx" "${cxx_std}"
    "--output" "${lock_file_path}"
)

if [ "${compiler}" != "default" ]; then
    cmd+=("--compiler-binary" "${compiler}")
fi

if [ -n "${flavor}" ]; then
    cmd+=("--flavor" "${flavor}")
fi

"${cmd[@]}"

checkpoint "Lock file prepared"

end_section



start_section "Create spack environment"
spack env create -d "${env_dir}" "${lock_file_path}" --with-view "$view_dir"
checkpoint "Spack environment created"
spack -e "${env_dir}" spec -l
checkpoint "Spack spec complete"
spack -e "${env_dir}" find
checkpoint "Spack find complete"
end_section

start_section "Install spack packages"
# Retry to absorb transient GHCR fetch failures (rate limits, network hiccups)
# when pulling binaries from the buildcache. These usually surface as spack
# reporting "no binary available" rather than as a network error, so that
# signature is retried too (see TRANSIENT_ERROR_PATTERNS). Install is
# idempotent: already-installed specs are skipped on subsequent attempts, so
# each retry only re-attempts what is still missing and is therefore cheap.
retry_transient spack -e "${env_dir}" install --fail-fast --use-buildcache only --concurrent-packages 10
checkpoint "Spack install complete"
end_section

start_section "Patch up Geant4 data directory"
if [ "${full_install:-false}" == "true" ]; then
  if ! which uv &> /dev/null ; then
    echo "uv not found, installing uv"
    curl -LsSf https://astral.sh/uv/install.sh | sh
    UV_EXE="/root/.local/bin/uv"
    checkpoint "uv installation complete"
  else
    UV_EXE=$(which uv)
  fi
  $UV_EXE run "$SCRIPT_DIR/download_geant4_datasets.py" -j8 --config "${view_dir}/bin/geant4-config"
  checkpoint "Geant4 datasets download complete"
fi
geant4_dir=$(spack -e "${env_dir}" location -i geant4)
# Prepare the folder for G4 data, and symlink it to where G4 will look for it.
# `data` itself, not just its parent: geant4.sh does `cd .../share/Geant4/data`
# and only the full_install path above creates it. Without it the ln below has
# no existing directory to resolve onto and writes a symlink at the target's own
# path -- a self-reference that makes every later cd fail with ELOOP.
mkdir -p "${geant4_dir}/share/Geant4/data"
[ -e "${view_dir}/share/Geant4/data" ] ||
  ln -s "${geant4_dir}/share/Geant4/data" "${view_dir}/share/Geant4/data"
end_section

start_section "Prepare python environment"
"${view_dir}/bin/python3" -m venv --system-site-packages "$venv_dir"
# NOTE: pip, not uv, on purpose. The venv is deliberately --system-site-packages
# so that the packages the spack view already provides (numpy and everything
# built against it) are reused rather than replaced. pip honours that and skips
# them; uv ignores system site-packages entirely and installs its own PyPI wheel
# over the top, which silently swaps out the spack-built stack.
retry_transient "${venv_dir}/bin/python3" -m pip install pyyaml jinja2
if [ "${full_install:-false}" == "true" ]; then
  retry_transient "${venv_dir}/bin/python3" -m pip install -r "${SCRIPT_DIR}/../../Python/Examples/tests/requirements.txt"
  retry_transient "${venv_dir}/bin/python3" -m pip install histcmp==0.10.0 matplotlib
  retry_transient "${venv_dir}/bin/python3" -m pip install pytest-md-report
fi
checkpoint "Python environment prepared"
end_section

start_section "Set environment variables"
set_env PATH "${venv_dir}/bin:${view_dir}/bin/:${PATH}"
set_env LD_LIBRARY_PATH "${venv_dir}/lib:${view_dir}/lib:${view_dir}/lib/root"
set_env DYLD_LIBRARY_PATH "${venv_dir}/lib:${view_dir}/lib:${view_dir}/lib/root"
set_env CMAKE_PREFIX_PATH "${venv_dir}:${view_dir}"
set_env ROOT_SETUP_SCRIPT "${view_dir}/bin/thisroot.sh"
set_env ROOT_INCLUDE_PATH "${view_dir}/include"
# cleanup setup-python mess
set_env PKG_CONFIG_PATH ""
set_env pythonLocation ""
set_env Python_ROOT_DIR ""
set_env Python2_ROOT_DIR ""
set_env Python3_ROOT_DIR ""
end_section

checkpoint "Setup script complete"
