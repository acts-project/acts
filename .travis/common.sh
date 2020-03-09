#!/bin/bash
set -e

echo "BEGIN: $(date)"

source CI/setup_lcg${LCG}.sh || true

mkdir build
cd build

cmake .. -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${COMMON_BUILD_OPTIONS}

cmake --build . -- -j$(nproc)

cmake --build . -- test
cmake --build . -- integrationtests
