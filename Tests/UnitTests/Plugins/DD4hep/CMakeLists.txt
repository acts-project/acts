
add_library(
    ActsTestsDD4hepFactories SHARED
    PrimitivesFactory.cpp)

target_include_directories(
    ActsTestsDD4hepFactories
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
    ActsTestsDD4hepFactories PUBLIC ActsPluginDD4hep DD4hep::DDCore)

set(ENV{LD_LIBRARY_PATH} ${CMAKE_INSTALL_PREFIX}/bin:$ENV{LD_LIBRARY_PATH})

set(unittest_extra_libraries ActsTestsDD4hepFactories ActsPluginDD4hep)
add_unittest(DD4hepDetectorElement DD4hepDetectorElementTests.cpp)
