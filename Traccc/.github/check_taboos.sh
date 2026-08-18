#!/bin/bash

if ! [ -x "$(command -v rg)" ]; then
   GREP="grep -R"
else
   GREP="rg -j 1 --no-heading -N"
fi

INPUT="$@"

INPUT_EX_THIS_FILE=${INPUT[@]/".github/check_taboos.sh"}

${GREP} "barcode" ${INPUT_EX_THIS_FILE[@]} ; test $? -eq 1 || exit 1

UNIT_CONSTANT_EXCUDE_FILE="core/include/traccc/definitions/common.hpp"
INPUT_EX_UNIT_CONSTANT=${INPUT_EX_THIS_FILE[@]/$UNIT_CONSTANT_EXCUDE_FILE}

${GREP} "detray::unit" ${INPUT_EX_UNIT_CONSTANT[@]} ; test $? -eq 1 || exit 1
${GREP} "detray::constant" ${INPUT_EX_UNIT_CONSTANT[@]} ; test $? -eq 1 || exit 1

TRACK_PARAMS_EXCLUDE_FILE="core/include/traccc/edm/track_parameters.hpp"
INPUT_EX_TRACK_PARAMS=${INPUT_EX_THIS_FILE[@]/$TRACK_PARAMS_EXCLUDE_FILE}

${GREP} "detray::free_track_parameters" ${INPUT_EX_TRACK_PARAMS[@]} ; test $? -eq 1 || exit 1
${GREP} "detray::bound_track_parameters" ${INPUT_EX_TRACK_PARAMS[@]} ; test $? -eq 1 || exit 1
${GREP} "detray::bound_vector" ${INPUT_EX_TRACK_PARAMS[@]} ; test $? -eq 1 || exit 1
${GREP} "detray::bound_matrix" ${INPUT_EX_TRACK_PARAMS[@]} ; test $? -eq 1 || exit 1
