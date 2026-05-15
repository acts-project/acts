#!/bin/bash

set -euo pipefail

function check_includes {
    file=$1
    out=$(printf "#include <${file:13}>\nint main() { return 0; }" | clang++ -std=c++20 -O0 -c -I "Core/include" -I "/usr/include/eigen3" -x c++ - 2>&1)
    if [[ "$?" -ne "0" ]]; then
        echo "$(date +%H:%M:%S)    Failed: $file"
        echo "$out"
        echo "------------------------------------"
    else
      # echo "$(date +%H:%M:%S)    Success: $file"
      :
    fi
}

export -f check_includes

find Core/include/ -type f -name "*.hpp" -o -name "*.ipp"| parallel check_includes
