#! /bin/bash

# bash3 workaround for dictionary

lcg_map_96=96
lcg_map_95=95apython3

lcg=$1
os=$2

i="lcg_map_${lcg}"
lcg_full="${!i}"

echo "LCG: $lcg -> $lcg_full"
echo "OS: $os"

set -e

source /opt/lcg/views/LCG_${lcg_full}/x86_64-centos7-gcc8-opt/setup.sh

mkdir build && cd build

cmake ..
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=17 \
  -DACTS_BUILD_UNITTESTS=ON \
  -DCMAKE_CXX_FLAGS="-Werror -fdiagnostics-color=always" \
  -DACTS_BUILD_DIGITIZATION_PLUGIN=on \
  -DACTS_BUILD_IDENTIFICATION_PLUGIN=on \
  -DACTS_BUILD_JSON_PLUGIN=on \
  -DACTS_BUILD_BENCHMARKS=on \
  -DACTS_BUILD_FATRAS=on \
  -DACTS_BUILD_EXAMPLES=on \
  -DACTS_BUILD_UNITTESTS=on \
  -DACTS_BUILD_LEGACY=on \
  -DACTS_BUILD_DD4HEP_PLUGIN=on \
  -DACTS_BUILD_TGEO_PLUGIN=on \
  -DACTS_BUILD_INTEGRATIONTESTS=on \

cmake --build . -- -j$(nproc)
cmake --build . -- test
cmake --build . -- integrationtests
