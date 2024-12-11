#!/bin/bash

set -u
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source ${SCRIPT_DIR}/detect_os.sh

packages_file=${GITHUB_WORKSPACE}/spack/etc/spack/packages.yaml

if [ "$os" == "ubuntu" ]; then
  sudo apt-get update
  sudo apt-get install -y libgl1-mesa-dev
cat <<EOF > "$packages_file"
packages:
  opengl:
    buildable: false
    externals:
    - prefix: /usr/
      spec: opengl@4.5
EOF
cat "$packages_file"
elif [ "$os" == "almalinux" ]; then
  dnf install -y mesa-libGLU
cat <<EOF > "$packages_file"
packages:
  opengl:
    buildable: false
    externals:
    - prefix: /usr/
      spec: opengl@4.6
EOF
cat "$packages_file"
fi
