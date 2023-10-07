#!/bin/bash
build_dir=$1

export NINJA_STATUS="[ninja][%f/%t] "

pushd $build_dir
ninja | grep -v '\[ninja\]'
popd
