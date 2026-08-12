# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include(traccc-functions)

# Specific flags for the NVIDIA backend of the HIP compiler.
if("${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia")
    traccc_add_flag( CMAKE_HIP_FLAGS "--expt-relaxed-constexpr" )
    traccc_add_flag( CMAKE_HIP_FLAGS "--use_fast_math" )
endif()

# Warnings and warnings-as-errors come from the shared ActsWarningFlags.cmake
# module. CMake picks the right flag for the platform in use (`-Werror` on amd,
# `-Werror all-warnings` on nvidia).
acts_apply_warning_flags(PROFILE traccc LANGUAGES HIP)
