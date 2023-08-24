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

# set test library dependencies
# set(unittest_extra_libraries ActsTestsDD4hepFactories ActsPluginDD4hep)
# add_unittest(DD4hepContainerStructure DD4hepContainerStructureTests.cpp)
# add_unittest(DD4hepDetectorElement DD4hepDetectorElementTests.cpp)
# add_unittest(DD4hepCylinderLayerStructure DD4hepCylinderLayerStructureTests.cpp)
# add_unittest(DD4hepDiscLayerStructure DD4hepDiscLayerStructureTests.cpp)
# add_unittest(DD4hepVolumeStructure DD4hepVolumeStructureTests.cpp)
