#!/bin/bash

set -e

output_dir=$1
shift
build_dir=$1
shift

mkdir -p $output_dir
output_dir=$(realpath $output_dir)

pushd $build_dir
NINJA_STATUS="[ninja] [%f/%t] " ninja $@ | tee $output_dir/ninja.log
popd

# grep fails if it does not find anything
set +e
cat $output_dir/ninja.log | grep -v '\[ninja\]' > $output_dir/clang-tidy.log
