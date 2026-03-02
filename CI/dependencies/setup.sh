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

# Parse command line arguments
while getopts "c:t:d:e:fh" opt; do
  case ${opt} in
    c )
      compiler=$OPTARG
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

echo "Ensure repo is synced with version ${_spack_repo_version}"

git config --global --add safe.directory "${_spack_repo_directory}"
spack repo update builtin --tag "${_spack_repo_version}"
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
    "--output" "${lock_file_path}"
)

if [ "${compiler}" != "default" ]; then
    cmd+=("--compiler-binary" "${compiler}")
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
spack -e "${env_dir}" install --fail-fast --use-buildcache only --concurrent-packages 10
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
# Prepare the folder for G4 data, and symlink it to where G4 will look for it
mkdir -p "${geant4_dir}/share/Geant4"
ln -s "${geant4_dir}/share/Geant4/data" "${view_dir}/share/Geant4/data"
end_section

start_section "Prepare python environment"
"${view_dir}/bin/python3" -m venv --system-site-packages "$venv_dir"
"${venv_dir}/bin/python3" -m pip install pyyaml jinja2
if [ "${full_install:-false}" == "true" ]; then
  "${venv_dir}/bin/python3" -m pip install -r "${SCRIPT_DIR}/../../Python/Examples/tests/requirements.txt"
  "${venv_dir}/bin/python3" -m pip install histcmp==0.8.2 matplotlib
  "${venv_dir}/bin/python3" -m pip install pytest-md-report
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
