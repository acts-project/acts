#!/bin/bash
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# Script generating TGZ/MD5 files that could then be uploaded to the ACTS web
# service.
#

# Stop on errors.
set -e
set -o pipefail

# Function printing the usage information for the script.
usage() {
   echo "Script generating data TGZ/MD5 files"
   echo ""
   echo "Usage: traccc_data_package_files.sh [options]"
   echo ""
   echo "Options:"
   echo "  -o <outputName>      Set the name of the output file(s)"
   echo "  -i <inputDirectory>  Additional input directory to pick up"
   echo "  -d <dataDirectory>   Main data directory"
   echo "  -c <cmakeExecutable> CMake executable to use in the script"
   echo ""
}

# Default script arguments.
TRACCC_DATA_NAME=${TRACCC_DATA_NAME:-"traccc-data-v11"}
TRACCC_DATA_DIRECTORY_NAMES=("cca_test" "detray_simulation" "geometries" "odd"
   "single_module" "tml_detector" "tml_full" "tml_pixel_barrel" "tml_pixels"
   "two_modules")
TRACCC_DATA_DIRECTORY=${TRACCC_DATA_DIRECTORY:-$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)}
TRACCC_CMAKE_EXECUTABLE=${TRACCC_CMAKE_EXECUTABLE:-cmake}

# Parse the command line argument(s).
while getopts ":o:i:d:ch" opt; do
   case $opt in
      o)
         TRACCC_DATA_NAME=$OPTARG
         ;;
      i)
         TRACCC_DATA_DIRECTORY_NAMES+=($OPTARG)
         ;;
      d)
         TRACCC_DATA_DIRECTORY=$OPTARG
         ;;
      c)
         TRACCC_CMAKE_EXECUTABLE=$OPTARG
         ;;
      h)
         usage
         exit 0
         ;;
      :)
         echo "Argument -$OPTARG requires a parameter!"
         usage
         exit 1
         ;;
      ?)
         echo "Unknown argument: -$OPTARG"
         usage
         exit 1
         ;;
   esac
done

# Go into the source directory.
cd "${TRACCC_DATA_DIRECTORY}"

# Compress the directories.
"${TRACCC_CMAKE_EXECUTABLE}" -E tar czf "${TRACCC_DATA_NAME}.tar.gz" \
   ${TRACCC_DATA_DIRECTORY_NAMES[@]}

# Generate an MD5 file.
"${TRACCC_CMAKE_EXECUTABLE}" -E md5sum "${TRACCC_DATA_NAME}.tar.gz" > \
   "${TRACCC_DATA_NAME}.md5"

# Leave the user with a message.
"${TRACCC_CMAKE_EXECUTABLE}" -E echo \
   "Generated files '${TRACCC_DATA_NAME}.tar.gz' and '${TRACCC_DATA_NAME}.md5'"
