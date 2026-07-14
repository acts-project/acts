#!/bin/bash
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# Script downloading the traccc data file(s) through HTTPS, and unpacking them.
#

# Stop on errors.
set -e
set -o pipefail

# Function printing the usage information for the script.
usage() {
   echo "Script downloading/unpacking data TGZ/MD5 files"
   echo ""
   echo "Usage: traccc_data_get_files.sh [options]"
   echo ""
   echo "Options:"
   echo "  -f <filename>        Name of the data file, without its extension"
   echo "  -d <webDirectory>    Directory holding the data and MD5 files"
   echo "  -o <dataDirectory>   Main data directory"
   echo "  -c <cmakeExecutable> CMake executable to use in the script"
   echo "  -w <curlExecutable>  CUrl executable to use in the script"
   echo ""
}

# Default script arguments.
TRACCC_DATA_NAME=${TRACCC_DATA_NAME:-"traccc-data-v10"}
TRACCC_WEB_DIRECTORY=${TRACCC_WEB_DIRECTORY:-"https://acts.web.cern.ch/traccc/data"}
TRACCC_DATA_DIRECTORY=${TRACCC_DATA_DIRECTORY:-$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)}
TRACCC_CMAKE_EXECUTABLE=${TRACCC_CMAKE_EXECUTABLE:-cmake}
TRACCC_CURL_EXECUTABLE=${TRACCC_CURL_EXECUTABLE:-curl}

# Parse the command line argument(s).
while getopts ":f:d:o:c:wh" opt; do
   case $opt in
      f)
         TRACCC_DATA_NAME=$OPTARG
         ;;
      d)
         TRACCC_WEB_DIRECTORY=$OPTARG
         ;;
      o)
         TRACCC_DATA_DIRECTORY=$OPTARG
         ;;
      c)
         TRACCC_CMAKE_EXECUTABLE=$OPTARG
         ;;
      w)
         TRACCC_CURL_EXECUTABLE=$OPTARG
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

# Go into the target directory.
cd "${TRACCC_DATA_DIRECTORY}"

# Download the TGZ and MD5 files.
"${TRACCC_CURL_EXECUTABLE}" --retry 5 --retry-delay 10 \
   --output "${TRACCC_DATA_NAME}.tar.gz"               \
   "${TRACCC_WEB_DIRECTORY}/${TRACCC_DATA_NAME}.tar.gz"
"${TRACCC_CURL_EXECUTABLE}" --retry 5 --retry-delay 10 \
   --output "${TRACCC_DATA_NAME}.md5"                  \
   "${TRACCC_WEB_DIRECTORY}/${TRACCC_DATA_NAME}.md5"

# Verify that the download succeeded.
"${TRACCC_CMAKE_EXECUTABLE}" -E md5sum "${TRACCC_DATA_NAME}.tar.gz" > \
   "${TRACCC_DATA_NAME}.md5-test"
"${TRACCC_CMAKE_EXECUTABLE}" -E compare_files "${TRACCC_DATA_NAME}.md5" \
   "${TRACCC_DATA_NAME}.md5-test"

# Extract the data files.
"${TRACCC_CMAKE_EXECUTABLE}" -E tar xf "${TRACCC_DATA_NAME}.tar.gz"

# Clean up.
"${TRACCC_CMAKE_EXECUTABLE}" -E remove "${TRACCC_DATA_NAME}.tar.gz" \
   "${TRACCC_DATA_NAME}.md5" "${TRACCC_DATA_NAME}.md5-test"

# Leave the user with a message.
"${TRACCC_CMAKE_EXECUTABLE}" -E echo \
   "Files from ${TRACCC_DATA_NAME} unpacked under '${TRACCC_DATA_DIRECTORY}'"
