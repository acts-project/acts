#!/bin/bash

set -e

build_dir=$1
shift

export NINJA_STATUS="[ninja][%f/%t] "

pushd $build_dir
ninja $@ | grep -v '\[ninja\]'
popd
