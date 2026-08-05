#!/bin/bash

if ! [ -x "$(command -v rg)" ]; then
   GREP="grep -R"
else
   GREP="rg -j 1 --no-heading -N"
fi

${GREP} "#include \"vecmem" "$@" ; test $? -eq 1 || exit 1
${GREP} "#include \"detray" "$@" ; test $? -eq 1 || exit 1
${GREP} "#include \"Acts" "$@" ; test $? -eq 1 || exit 1
${GREP} "#include \"covfie" "$@" ; test $? -eq 1 || exit 1
${GREP} "#include \"algebra" "$@" ; test $? -eq 1 || exit 1
