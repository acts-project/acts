# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

message(STATUS "Building the SMatrix algebra plugin")

find_package(ROOT COMPONENTS Core MathCore Smatrix REQUIRED)

# Set up the library.
detray_add_library( detray_algebra_smatrix algebra_smatrix
   "include/algebra/smatrix.hpp"
   "include/algebra/impl/*.hpp"
   "include/algebra/impl/detail/*.hpp"
)

target_link_libraries(
    detray_algebra_smatrix
    INTERFACE ROOT::Core ROOT::MathCore ROOT::Smatrix detray::core
)

target_compile_definitions(
    detray_algebra_smatrix
    INTERFACE
        DETRAY_CUSTOM_SCALARTYPE=${DETRAY_CUSTOM_SCALARTYPE}
        DETRAY_ALGEBRA_SMATRIX
)

detray_test_public_headers( detray_algebra_smatrix
   "algebra/smatrix.hpp"
)
