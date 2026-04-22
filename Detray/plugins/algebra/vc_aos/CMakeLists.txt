# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the Vc AoS algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_vc_aos algebra_vc_aos
   "include/algebra/vc_aos.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(detray_algebra_vc_aos INTERFACE Vc::Vc detray::core)

target_compile_definitions(
    detray_algebra_vc_aos
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_VC_AOS
)

detray_test_public_headers( detray_algebra_vc_aos
   "algebra/vc_aos.hpp"
)
