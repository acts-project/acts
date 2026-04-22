# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the std::array algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_array algebra_array
   "include/algebra/array.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(detray_algebra_array INTERFACE detray::core)

target_compile_definitions(
    detray_algebra_array
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_ARRAY
)

detray_test_public_headers( detray_algebra_array
   "algebra/array.hpp"
)
