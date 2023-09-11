#!/bin/bash

scriptdir="$( dirname -- "$BASH_SOURCE"; )";

input=`find "$1" -iname *.cpp -or -iname *.hpp -or -iname *.ipp -not -path "thirdparty/*"`

python "$scriptdir/check_end_of_file.py" --reject-multiple-newlines $input
