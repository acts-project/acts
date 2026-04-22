# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the Vc SoA algebra plugin")

# Set up the library.
detray_add_library( detray_algebra_vc_soa algebra_vc_soa
   "include/algebra/vc_soa.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/utils/*.hpp"
)

target_link_libraries(
    detray_algebra_vc_soa
    INTERFACE Vc::Vc detray::algebra_vc_aos
)

target_compile_definitions(
    detray_algebra_vc_soa
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_VC_SOA
)

detray_test_public_headers( detray_algebra_vc_soa
   "algebra/vc_soa.hpp"
)
