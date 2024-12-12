#!/bin/bash
set -e
set -u

source "$(spack location -r)/share/spack/setup-env.sh"

env=$1
shift
spack env activate "$env"

"$@"
