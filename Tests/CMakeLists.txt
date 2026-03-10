# helpers shared between all tests
add_subdirectory(CommonHelpers)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(
        "-Wno-c2y-extensions"
        ACTS_TESTS_HAS_WNO_C2Y_EXTENSIONS
    )
endif()

add_subdirectory_if(Benchmarks ACTS_BUILD_BENCHMARKS)
add_subdirectory_if(IntegrationTests ACTS_BUILD_INTEGRATIONTESTS)
add_subdirectory_if(UnitTests ACTS_BUILD_UNITTESTS)
