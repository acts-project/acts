#!/bin/bash
set -e
set -u

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

export SPACK_COLOR=always

function set_env {
  key="$1"
  value="$2"

  echo "=> ${key}=${value}"

  if [ -n "${GITHUB_ACTIONS:-}" ]; then
    echo "${key}=${value}" >> "$GITHUB_ENV"
  else
    export "${key}"="${value}"
  fi
}

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

tag=${1:-${DEPENDENCY_TAG}}
if [ -n "${GITHUB_ACTIONS:-}" ]; then
    destination="${GITHUB_WORKSPACE}/dependencies"
elif [ -n "${GITLAB_CI:-}" ];then
    destination="${CI_PROJECT_DIR}/dependencies"
else
    destination=${2}
fi

echo "Install tag: $tag"
echo "Install destination: $destination"

mkdir -p ${destination}


start_section "Install spack if not already installed"
if ! command -v spack &> /dev/null; then
    git clone -c feature.manyFiles=true https://github.com/spack/spack.git
    pushd spack > /dev/null
    git config user.name github-actions[bot]
    git config user.email 41898282+github-actions[bot]@users.noreply.github.com
    # Apply patch for spack improvements
    curl https://patch-diff.githubusercontent.com/raw/spack/spack/pull/47370.patch | git am
    source "$(pwd)/share/spack/setup-env.sh"
    popd > /dev/null
fi
end_section

if [ -n "${CI:-}" ]; then
start_section "Add buildcache mirror"
mirror_name="acts-spack-buildcache"
mirror_url="oci://ghcr.io/acts-project/spack-buildcache"
spack mirror add ${mirror_name} ${mirror_url} --unsigned
end_section

start_section "Locate OpenGL"
"${SCRIPT_DIR}/opengl.sh"
end_section
fi

start_section "Get spack lock file name"
arch=$(spack arch --family)
lock_file="spack-${arch}.lock"
end_section


url="https://github.com/acts-project/ci-dependencies/releases/download/${tag}/${lock_file}"
echo "URL: $url"

env_dir="${destination}/env"
view_dir="${destination}/view"
mkdir -p ${env_dir}
lock_file_path="${destination}/spack.lock"
start_section "Get spack lock files from tag"
curl -fL -o $lock_file_path $url
end_section



start_section "Create spack environment"
time spack env create -d "${env_dir}" "${lock_file_path}" --with-view "$view_dir"
time spack -e "${env_dir}" spec -l
time spack -e "${env_dir}" find
end_section

start_section "Install spack packages"
if [ "$(uname)" = "Darwin" ]; then
  NCPUS=$(sysctl -n hw.ncpu)
else
  NCPUS=$(nproc)
fi
time "${SCRIPT_DIR}"/parallel.sh "$NCPUS" spack -e "${env_dir}" install --use-buildcache only \
  | tee install.log \
  | grep -v "^Waiting\|^\[+\]"
end_section

start_section "Patch up Geant4 data directory"
# ${SCRIPT_DIR}/with_spack_env.sh ${env_dir} geant4-config --install-datasets
geant4_dir=$(spack -e "${env_dir}" location -i geant4)
# Prepare the folder for G4 data, and symlink it to where G4 will look for it
mkdir -p "${geant4_dir}"/share/Geant4
ln -s "${geant4_dir}"/share/Geant4/data ${view_dir}/share/Geant4/data
end_section

start_section "Set environment variables"
set_env CMAKE_PREFIX_PATH "${view_dir}"
set_env LD_LIBRARY_PATH "${view_dir}/lib"
set_env ROOT_INCLUDE_PATH "${view_dir}/include"
# Geant4 puts CLHEP in a subdirectory
set_env ROOT_INCLUDE_PATH "${view_dir}/include/Geant4"
end_section

start_section "Prepare python environment"
ls -al
venv_dir="${view_dir}/venv"
"${view_dir}"/bin/python3 -m venv "$venv_dir"

if [ -n "${GITHUB_ACTIONS:-}" ]; then
  echo "${venv_dir}/bin" >> "$GITHUB_PATH"
  echo "${view_dir}/bin" >> "$GITHUB_PATH"
fi
set_env PATH "${venv_dir}/bin:${view_dir}/bin/:${PATH}"
end_section



# Pythia8 looks for settings in this directory
# set_env PYTHIA8DATA "${destination}/share/Pythia8/xmldoc"
