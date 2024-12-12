#!/bin/bash
set -e
set -u

source "${SPACK_ROOT}/share/spack/setup-env.sh"

env=$1
shift
spack env activate "$env"

"$@"
