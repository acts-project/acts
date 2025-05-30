#!/bin/bash

RET=0
ERRORS=0

FILES=$(find Core/include/ -type f -name "*.hpp" -o -name "*.ipp")
N_FILES=$(echo "$FILES" | wc -l)
echo "Check $N_FILES files"

ITER=0

for file in $FILES; do
    ITER=$((ITER+1))
    echo "$(date +%H:%M:%S)    $((100*ITER/N_FILES))%   check $file"
    out=$(printf "#include <${file:13}>\nint main() { return 0; }" | clang++ -std=c++20 -O0 -c -I "Core/include" -I "/usr/include/eigen3" -x c++ - 2>&1)
    if [[ "$?" -ne "0" ]]; then
        echo "------------------------------------"
        echo "$out"
        echo "------------------------------------"
        RET=1
        ERRORS=$((ERRORS+1))
    fi
done

echo "Total errors: $ERRORS"
exit $RET
