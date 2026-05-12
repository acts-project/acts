#Detray library, part of the ACTS project(R& D line)
#
#(c) 2023-2025 CERN for the benefit of the ACTS project
#
#Mozilla Public License Version 2.0

message(STATUS "Building detray host command line tools")

include(CMakeFindDependencyMacro)

find_dependency(Boost COMPONENTS program_options)

# Generate test detectors
detray_add_executable(generate_toy_detector
                      "generate_toy_detector.cpp"
                      LINK_LIBRARIES Boost::program_options detray::tools
                      detray::io detray::test_common detray::test_utils detray::core_array
)

detray_add_executable(generate_wire_chamber
                      "generate_wire_chamber.cpp"
                      LINK_LIBRARIES Boost::program_options detray::tools
                      detray::io detray::test_common detray::test_utils detray::core_array
)

detray_add_executable(generate_telescope_detector
                      "generate_telescope_detector.cpp"
                      LINK_LIBRARIES Boost::program_options detray::tools
                      detray::io detray::test_common detray::test_utils detray::core_array
)

if(DETRAY_SVG_DISPLAY)
    # Build the visualization executable.
    detray_add_executable(detector_display
                        "detector_display.cpp"
                        LINK_LIBRARIES Boost::program_options detray::core_array detray::io detray::tools
                        detray::svgtools
    )
endif()

if(DETRAY_BUILD_TESTING)
    # Build the navigation validation executable.
    detray_add_executable(navigation_validation
                        "navigation_validation.cpp"
                        LINK_LIBRARIES GTest::gtest GTest::gtest_main
                        Boost::program_options detray::core_array detray::test_cpu detray::tools
    )

    # Build the material validation executable.
    detray_add_executable(material_validation
                        "material_validation.cpp"
                        LINK_LIBRARIES GTest::gtest GTest::gtest_main
                        Boost::program_options detray::core_array detray::test_cpu detray::tools
    )
endif()

if(DETRAY_BUILD_BENCHMARKS)
    # Look for openMP, which is used for the CPU propagation benchmark
    find_package(OpenMP)

    # Build the propagation benchmark executable.
    macro(detray_add_propagation_benchmark algebra type)
        detray_add_executable(propagation_${type}_cpu_${algebra}
                            "propagation_${type}.cpp"
                            LINK_LIBRARIES Boost::program_options benchmark::benchmark benchmark::benchmark_main vecmem::core detray::benchmark_cpu detray::core_${algebra} detray::tools detray::detectors
        )

        if(OpenMP_CXX_FOUND)
            target_link_libraries(
                detray_propagation_${type}_cpu_${algebra}
                PRIVATE OpenMP::OpenMP_CXX
            )
        endif()
    endmacro()

    # Build the array benchmark.
    detray_add_propagation_benchmark( array benchmark )
    detray_add_propagation_benchmark( array scaling )

    # Build the Eigen benchmark executable.
    if(DETRAY_EIGEN_PLUGIN)
        detray_add_propagation_benchmark( eigen benchmark )
        detray_add_propagation_benchmark( eigen scaling )
    endif()

    # Build the SMatrix benchmark executable.
    if(DETRAY_SMATRIX_PLUGIN)
        detray_add_propagation_benchmark( smatrix benchmark )
        detray_add_propagation_benchmark( smatrix scaling )
    endif()

    # Build the Vc benchmark executable.
    if(DETRAY_VC_AOS_PLUGIN)
        detray_add_propagation_benchmark( vc_aos benchmark )
        detray_add_propagation_benchmark( vc_aos scaling )
    endif()
endif()
