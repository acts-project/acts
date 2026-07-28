# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include(traccc-functions)

# Only tweak the flags for the Intel compiler.
if(
    NOT (
        ("${CMAKE_SYCL_COMPILER_ID}" STREQUAL "IntelLLVM")
        OR ("${CMAKE_SYCL_COMPILER_ID}" MATCHES "Clang")
    )
)
    return()
endif()

# All SYCL warning flags come from the shared ActsWarningFlags.cmake module.
# Note that CMake has no COMPILE_WARNING_AS_ERROR implementation for SYCL (it is
# a vecmem-provided pseudo-language), so the module emits -Werror by hand.
acts_apply_warning_flags(PROFILE traccc LANGUAGES SYCL)
