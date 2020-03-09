#!/bin/bash
set -e

echo "BEGIN: $(date)"

source CI/setup_lcg96.sh || true

mkdir build
cd build


cmake .. -GNinja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS="-Werror -fdiagnostics-color=always" \
  -DACTS_BUILD_DIGITIZATION_PLUGIN=on \
  -DACTS_BUILD_IDENTIFICATION_PLUGIN=on \
  -DACTS_BUILD_JSON_PLUGIN=on \
  -DACTS_BUILD_BENCHMARKS=on \
  -DACTS_BUILD_FATRAS=on \
  -DACTS_BUILD_EXAMPLES=on \
  -DACTS_BUILD_INTEGRATIONTESTS=on \
  -DACTS_BUILD_UNITTESTS=on

cmake --build . -- -j$(nproc)

export CTEST_OUTPUT_ON_FAILURE=1

cmake --build . -- test
