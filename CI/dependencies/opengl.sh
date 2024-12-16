#!/bin/bash

set -u
set -e

packages_file=$(spack location -r)/etc/spack/packages.yaml
echo "Packages file: $packages_file"
stat "$packages_file" || true

if ! command -v sudo &> /dev/null
then
    SUDO=""
else
    SUDO="sudo"
fi

os=$(spack arch --family)

echo "OS: $os"

if [[ "$os" == *ubuntu* ]]; then
  ${SUDO} apt-get update
  ${SUDO} apt-get install -y libgl1-mesa-dev

if [[ "$os" == *ubuntu24* ]]; then
  version="4.6"
elif [[ "$os" == *ubuntu20* ]]; then
  version="4.5"
else
  echo "Unknown OS version, default OpenGL version"
  version="4.5"
fi

cat <<EOF > "$packages_file"
packages:
  opengl:
    buildable: false
    externals:
    - prefix: /usr/
      spec: opengl@${version}
EOF
cat "$packages_file"
elif [[ "$os" == *almalinux* ]]; then
  ${SUDO} dnf install -y mesa-libGLU
cat <<EOF > "$packages_file"
packages:
  opengl:
    buildable: false
    externals:
    - prefix: /usr/
      spec: opengl@4.6
EOF
cat "$packages_file"
else [[ "$os" == *darwin* ]]
  echo "Nothing to do on Darwin"
fi
