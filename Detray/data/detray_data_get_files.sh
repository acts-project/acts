#!/bin/bash
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0
#
# Script downloading the detray data file(s) through HTTPS, and unpacking them.
#

# Stop on errors.
set -e
set -o pipefail

# Function printing the usage information for the script.
usage() {
   echo "Script downloading/unpacking data TGZ/MD5 files"
   echo ""
   echo "Usage: detray_data_get_files.sh [options]"
   echo ""
   echo "Options:"
   echo "  -f <filename>        Name of the data file, without its extension"
   echo "  -d <webDirectory>    Directory holding the data and MD5 files"
   echo "  -o <dataDirectory>   Main data directory"
   echo "  -c <cmakeExecutable> CMake executable to use in the script"
   echo "  -w <curlExecutable>  Curl executable to use in the script"
   echo ""
   return 0
}

# Default script arguments.
DETRAY_DATA_NAME=${DETRAY_DATA_NAME:-"detray-data-v2"}
DETRAY_WEB_DIRECTORY=${DETRAY_WEB_DIRECTORY:-"https://acts.web.cern.ch/traccc/data"}
DETRAY_DATA_DIRECTORY=${DETRAY_DATA_DIRECTORY:-$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)}
DETRAY_CMAKE_EXECUTABLE=${DETRAY_CMAKE_EXECUTABLE:-cmake}
DETRAY_CURL_EXECUTABLE=${DETRAY_CURL_EXECUTABLE:-curl}

# Parse the command line argument(s).
while getopts ":f:d:o:c:wh" opt; do
   case $opt in
      f)
         DETRAY_DATA_NAME=$OPTARG
         ;;
      d)
         DETRAY_WEB_DIRECTORY=$OPTARG
         ;;
      o)
         DETRAY_DATA_DIRECTORY=$OPTARG
         ;;
      c)
         DETRAY_CMAKE_EXECUTABLE=$OPTARG
         ;;
      w)
         DETRAY_CURL_EXECUTABLE=$OPTARG
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
      *)
         echo "Unknown argument: -$OPTARG"
         usage
         exit 1
         ;;
   esac
done

# Go into the target directory.
cd "${DETRAY_DATA_DIRECTORY}"

# Download the TGZ and MD5 files.
"${DETRAY_CURL_EXECUTABLE}" --retry 5 --retry-delay 10 \
   --output "${DETRAY_DATA_NAME}.tar.gz"               \
   "${DETRAY_WEB_DIRECTORY}/${DETRAY_DATA_NAME}.tar.gz"
"${DETRAY_CURL_EXECUTABLE}" --retry 5 --retry-delay 10 \
   --output "${DETRAY_DATA_NAME}.md5"                  \
   "${DETRAY_WEB_DIRECTORY}/${DETRAY_DATA_NAME}.md5"

# Verify that the download succeeded.
"${DETRAY_CMAKE_EXECUTABLE}" -E md5sum "${DETRAY_DATA_NAME}.tar.gz" > \
   "${DETRAY_DATA_NAME}.md5-test"
"${DETRAY_CMAKE_EXECUTABLE}" -E compare_files "${DETRAY_DATA_NAME}.md5" \
   "${DETRAY_DATA_NAME}.md5-test"

# Extract the data files.
"${DETRAY_CMAKE_EXECUTABLE}" -E tar xf "${DETRAY_DATA_NAME}.tar.gz"

# Clean up.
"${DETRAY_CMAKE_EXECUTABLE}" -E remove "${DETRAY_DATA_NAME}.tar.gz" \
   "${DETRAY_DATA_NAME}.md5" "${DETRAY_DATA_NAME}.md5-test"

# Leave the user with a message.
"${DETRAY_CMAKE_EXECUTABLE}" -E echo \
   "Files from ${DETRAY_DATA_NAME} unpacked under '${DETRAY_DATA_DIRECTORY}'"
