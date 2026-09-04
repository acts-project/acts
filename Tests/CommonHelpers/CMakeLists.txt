# get the path of the test data directory
file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../Data" acts_test_data_dir)

add_library(
    ActsTestsCommonHelpers
    SHARED
    src/DataDirectory.cpp
    src/DetectorElementStub.cpp
    src/CubicTrackingGeometry.cpp
    src/CylindricalTrackingGeometry.cpp
)

add_library(Acts::TestsCommonHelpers ALIAS ActsTestsCommonHelpers)

target_compile_definitions(
    ActsTestsCommonHelpers
    PRIVATE
        "ACTS_TEST_DATA_DIR=\"${acts_test_data_dir}\""
        BOOST_FILESYSTEM_NO_DEPRECATED
)
target_include_directories(
    ActsTestsCommonHelpers
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
if(ACTS_BUILD_EXAMPLES)
    target_link_libraries(
        ActsTestsCommonHelpers
        PUBLIC Acts::Core Acts::ExamplesFramework
    )
else()
    target_link_libraries(ActsTestsCommonHelpers PUBLIC ActsCore)
endif()

acts_compile_headers(TestsCommonHelpers GLOB "include/ActsTests/**/*.hpp")
