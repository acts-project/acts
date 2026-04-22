# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up all of the "host" tutorials.
if(DETRAY_BUILD_HOST)
    add_subdirectory(cpu)
endif()

# Set up all of the "device" tutorials.
if(DETRAY_BUILD_CUDA OR DETRAY_BUILD_SYCL)
    add_subdirectory(device)
endif()
