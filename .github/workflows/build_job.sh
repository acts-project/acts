#! /bin/bash

# bash3 workaround for dictionary

lcg_map_96=96
lcg_map_95=95apython3

os_map_cc7=centos7
os_map_slc6=slc6

lcg=$1
os=$2

i="lcg_map_${lcg}"
_lcg="${!i}"
i="os_map_${os}"
_os="${!i}"

echo "LCG: $lcg -> $_lcg"
echo "OS: $os -> $_os"


source /opt/lcg/views/LCG_${_lcg}/x86_64-${_os}-gcc8-opt/setup.sh

set -e

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
