# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# This is a simple wrapper around the ROOTConfig.cmake file that would come
# with a ROOT installation. It serves two purposes:
#  - Provides a nice printout about ROOT having been found, which
#    ROOTConfig.cmake doesn't do for some reason;
#  - Avoids printing scary-looking warnings from CMake when ROOT is
#    not found.

# Look for ROOTConfig.cmake.
find_package(ROOT CONFIG QUIET)

# Print a standard output about ROOT (not) having been found.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    ROOT
    FOUND_VAR ROOT_FOUND
    REQUIRED_VARS ROOT_INCLUDE_DIRS ROOT_LIBRARY_DIR ROOT_BINDIR
    VERSION_VAR ROOT_VERSION
)
