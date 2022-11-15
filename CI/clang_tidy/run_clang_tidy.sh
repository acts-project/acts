#!/bin/bash
build_dir=$1

pushd $build_dir
ninja -v
popd
