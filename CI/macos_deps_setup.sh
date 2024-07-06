#!/bin/bash
set -e
set -u

CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-}

function add_dep() {
  pkg=$1
  ver=$2
  if [ -z "$CMAKE_PREFIX_PATH" ]; then
    export CMAKE_PREFIX_PATH="$DEPENDENCY_DIR/$pkg/$ver"
  else
    export CMAKE_PREFIX_PATH="$DEPENDENCY_DIR/$pkg/$ver:$CMAKE_PREFIX_PATH"
  fi
}

add_dep python 3.12.2
add_dep tbb 2021.11.0
add_dep geant4 11.1.3
add_dep hepmc3 3.2.5
add_dep pythia8 311
add_dep nlohmann_json 3.11.2
add_dep root 6.30.06
add_dep podio 00-17-02
add_dep edm4hep 00-10-01
add_dep dd4hep 01-27
add_dep boost 1.84.0
add_dep eigen 3.4.0

echo $CMAKE_PREFIX_PATH
