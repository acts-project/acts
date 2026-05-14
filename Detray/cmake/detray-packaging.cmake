# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# CMake include(s).
include(CPack)

# Export the configuration of the project.
include(CMakePackageConfigHelpers)
set(CMAKE_INSTALL_CMAKEDIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/detray-${PROJECT_VERSION}"
)
install(
    EXPORT detray-exports
    NAMESPACE "detray::"
    FILE "detray-config-targets.cmake"
    DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/detray-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/detray-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
    PATH_VARS
        CMAKE_INSTALL_INCLUDEDIR
        CMAKE_INSTALL_LIBDIR
        CMAKE_INSTALL_CMAKEDIR
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/detray-config-version.cmake"
    COMPATIBILITY "AnyNewerVersion"
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/detray-config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/detray-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
)
