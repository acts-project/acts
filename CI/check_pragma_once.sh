#!/bin/bash

ec=0

for file in $(find Core Examples Tests Plugins -name "*.hpp"); do
    res=$(grep -e "^[[:space:]]*#pragma once" $file)
    if [[ "$res" != "#pragma once" ]]; then
        ec=1
        echo "'#pragma once' missing in '$file'"
    fi
done

exit $ec
