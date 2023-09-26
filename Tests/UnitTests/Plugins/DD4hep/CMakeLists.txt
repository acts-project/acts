add_library(
    ActsTestsDD4hepFactories SHARED
    DD4hepTestsHelper.cpp
    LayerFactory.cpp
    PrimitivesFactory.cpp)

target_include_directories(
    ActsTestsDD4hepFactories
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
    ActsTestsDD4hepFactories PUBLIC ActsPluginDD4hep DD4hep::DDCore)


set(DD4HEP_LIBRARY_PATH "${DD4hep_DIR}/lib")

# set test library dependencies
set(unittest_extra_libraries ActsTestsDD4hepFactories ActsPluginDD4hep)
add_unittest(DD4hepDetectorElement DD4hepDetectorElementTests.cpp)
if(APPLE)
    set_property(TEST DD4hepDetectorElement PROPERTY ENVIRONMENT "DYLD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH}")
else()
    set_property(TEST DD4hepDetectorElement PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${DD4HEP_LIBRARY_PATH}")
endif()