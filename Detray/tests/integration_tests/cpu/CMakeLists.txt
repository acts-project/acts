# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Setup integration tests that don't depend on algebra plugins
detray_add_integration_test(cpu
       "utils/ranges.cpp"
       LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core
)

# Run integration tests only after unit tests passed
set_tests_properties(
    detray_integration_test_cpu
    PROPERTIES DEPENDS "detray_unit_test_cpu"
)

# Macro setting up the CPU tests for a specific algebra plugin.
macro(detray_add_cpu_test algebra)
    # Build the test executable.
    detray_add_integration_test(cpu_${algebra}
       "builders/grid_builder.cpp"
       "builders/homogeneous_material_builder.cpp"
       "builders/material_map_builder.cpp"
       "builders/volume_builder.cpp"
       "material/material_interaction.cpp"
       "propagator/backward_propagation.cpp"
       "propagator/covariance_transport.cpp"
       "propagator/guided_navigator.cpp"
       "propagator/propagation_caching_navigator.cpp"
       "propagator/propagation_direct_navigator.cpp"
       "utils/type_registry.cpp"
       LINK_LIBRARIES GTest::gtest GTest::gtest_main detray::core_${algebra}
                      vecmem::core detray::test_framework detray::test_common detray::test_utils
    )

    # Run integration tests only after unit tests passed
    set_tests_properties(
        detray_integration_test_cpu_${algebra}
        PROPERTIES
            DEPENDS
                "detray_unit_test_cpu;detray_unit_test_covfie;detray_unit_test_cpu_${algebra}"
    )
endmacro()

# Build the array tests.
detray_add_cpu_test( array )

# Build the Eigen tests.
if(DETRAY_EIGEN_PLUGIN)
    detray_add_cpu_test( eigen )
endif()

# Build the Fastor tests.
if(DETRAY_FASTOR_PLUGIN)
    detray_add_cpu_test( fastor )
endif()

# Build the SMatrix tests.
if(DETRAY_SMATRIX_PLUGIN)
    detray_add_cpu_test( smatrix )
endif()

# Build the Vc AoS tests.
if(DETRAY_VC_AOS_PLUGIN)
    detray_add_cpu_test( vc_aos )
endif()

# Detector navigation tests
detray_add_integration_test(telescope_detector
                            "detectors/telescope_detector_navigation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            detray::test_cpu detray::io detray::test_framework detray::test_common
                            detray::core_array
)

set_tests_properties(
    detray_integration_test_telescope_detector
    PROPERTIES
        DEPENDS
            "detray_unit_test_cpu;detray_unit_test_cpu_array;detray_unit_test_svgtools"
)

detray_add_integration_test(toy_detector
                            "detectors/toy_detector_navigation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            detray::test_cpu detray::io detray::test_framework  detray::test_common
                            detray::core_array
)

set_tests_properties(
    detray_integration_test_toy_detector
    PROPERTIES
        DEPENDS
            "detray_unit_test_cpu;detray_unit_test_cpu_array;detray_unit_test_svgtools"
)

detray_add_integration_test(wire_chamber
                            "detectors/wire_chamber_navigation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            detray::test_cpu detray::io detray::test_framework detray::test_common
                            detray::core_array
)

set_tests_properties(
    detray_integration_test_wire_chamber
    PROPERTIES
        DEPENDS
            "detray_unit_test_cpu;detray_unit_test_cpu_array;detray_unit_test_svgtools"
)

include(CMakeFindDependencyMacro)

find_dependency(Boost COMPONENTS program_options)

detray_add_integration_test(jacobian_validation
                            "propagator/jacobian_validation.cpp"
                            LINK_LIBRARIES GTest::gtest GTest::gtest_main
                            Boost::program_options detray::core_array
                            detray::test_framework detray::test_common detray::test_utils detray::tools
)

set_tests_properties(
    detray_integration_test_jacobian_validation
    PROPERTIES
        DEPENDS
            "detray_unit_test_cpu;detray_unit_test_cpu_arraydetray_integration_test_telescope_detector"
)
