# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

if(DETRAY_BUILD_CUDA)
    add_subdirectory(cuda)
endif()

if(DETRAY_BUILD_HIP)
    add_subdirectory(hip)
endif()
