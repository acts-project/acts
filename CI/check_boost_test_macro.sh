#!/bin/bash

test_string="BOOST_TEST("
grep $test_string -n -r Tests --include "*.cpp" --include "*.hpp" --include "*.ipp"

status=$?

if [[ $status -eq 0 ]]; then
  echo "Found occurrences of '$test_string'"
  exit 1
else
  echo "Did not find occurrences of '$test_string'"
  exit 0
fi
