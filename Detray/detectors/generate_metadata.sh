# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

#!/bin/sh

# Print help message
print_help()
{
   echo "Syntax: generate_metadata.sh [-m|b|o|v|vv|h]"
   echo "Options:"
   echo "  m     Full path and filename of the python metadata script"
   echo "        in the \"detray/detectors/python\" directory."
   echo "  p     Where to find the detray python package (e.g. build/python)."
   echo "  o     Output directory."
   echo "  v     Verbose logging."
   echo "  vv    Debug level logging."
   echo "  h     Print this help."
   echo

   return
}

# Verbosity level of the generator script
log_lvl=0

# Parse options
while getopts ":hvp:vv:m:o:" opt; do
    case ${opt} in
        h)
            echo "Generate custom detray Detector Metadata"
            echo
            print_help
            exit 0
        ;;
        v)
            ((log_lvl++))
        ;;
        p)
            python_dir=${OPTARG}
        ;;
        m)
            metadata_generator=${OPTARG}
        ;;
        o)
            out_dir=${OPTARG}
        ;;
        \?)
            >&2 echo
            >&2 echo "ERROR: Invalid option ${opt}! Usage:"
            >&2 echo
            print_help
            exit 1
        ;;
        *)
            >&2 echo
            >&2 echo "ERROR: Unknown option ${opt}! Usage:"
            >&2 echo
            print_help
            exit 1
        ;;
    esac
done

if [ -z "${metadata_generator}" ];then
    >&2 echo
    >&2 echo "ERROR: No detector metadata script supplied! Usage:"
    >&2 echo
    print_help
    exit 1
fi

if [ -z "${python_dir}" ];then
    >&2 echo
    >&2 echo "ERROR: Path to the detray python package not supplied! Usage:"
    >&2 echo
    print_help
    exit 1
fi

# Check how to invoke the python interpreter
python_command="python3"
if command -v {$python_command} > /dev/null; then
    # Try again
    python_command="python"
    if command -v {$python_command} > /dev/null; then
        >&2 echo
        >&2 echo "ERROR: No python command found: exiting"
        >&2 echo
        exit 2
    fi
fi

generator_command="${python_command} ${metadata_generator} --no-format"

# Configure verbosity
if [ ${log_lvl} -eq 1 ]; then
    generator_command="${generator_command} -v"
elif [ ${log_lvl} -eq 2 ]; then
    generator_command="${generator_command} -vv"
fi

# Add the option for the custom output location
if [ ! -z "${out_dir}" ]; then
    generator_command="${generator_command} -o ${out_dir}"
fi

# Set up the detray python package
export PYTHONPATH="${python_dir}/python":$PYTHONPATH

# Run the generation of the requested metadata
eval ${generator_command}
