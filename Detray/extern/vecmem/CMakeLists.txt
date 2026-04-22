# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Silence FetchContent warnings with CMake >=3.24.
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Tell the user what's happening.
message(STATUS "Building VecMem as part of the Detray project")

# Declare where to get VecMem from.
set(DETRAY_VECMEM_SOURCE
    "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v${DETRAY_VECMEM_VERSION}.tar.gz;URL_MD5;${DETRAY_VECMEM_DIGEST}"
    CACHE STRING
    "Source for VecMem, when built as part of this project"
)
mark_as_advanced(DETRAY_VECMEM_SOURCE)

# Mark the import as a system library on modern CMake versions
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25.0)
    set(DETRAY_VECMEM_SOURCE_FULL "${DETRAY_VECMEM_SOURCE};SYSTEM")
else()
    set(DETRAY_VECMEM_SOURCE_FULL "${DETRAY_VECMEM_SOURCE}")
endif()
mark_as_advanced(DETRAY_VECMEM_SOURCE_FULL)

FetchContent_Declare(VecMem ${DETRAY_VECMEM_SOURCE_FULL})

# Options used in the build of VecMem.
set(VECMEM_BUILD_TESTING
    FALSE
    CACHE BOOL
    "Turn off the build of the VecMem unit tests"
)

# Get it into the current directory.
FetchContent_MakeAvailable(VecMem)
