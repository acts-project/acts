#!/bin/bash

set -e

output_dir=$1
shift
build_dir=$1
shift

mkdir -p $output_dir
output_dir=$(realpath $output_dir)

pushd $build_dir
print "why?"
NINJA_STATUS="[ninja][%f/%t] " ninja $@ | tee $output_dir/ninja.log
print "what?"
popd
print "hi!"
cat $output_dir/ninja.log | grep -v '\[ninja\]' > $output_dir/clang-tidy.log

print "hello?"
exit 0
