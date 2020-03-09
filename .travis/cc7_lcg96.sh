#!/bin/bash
set -e

echo "BEGIN: $(date)"

source CI/setup_lcg96.sh || true

mkdir build
cd build


cmake .. -DCMAKE_BUILD_TYPE=Release ${COMMON_BUILD_OPTIONS}

cmake --build . -- -j$(nproc)

export CTEST_OUTPUT_ON_FAILURE=1

cmake --build . -- test
