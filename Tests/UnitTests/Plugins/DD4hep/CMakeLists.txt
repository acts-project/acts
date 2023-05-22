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

# set test en    
set(unittest_extra_libraries ActsTestsDD4hepFactories ActsPluginDD4hep ActsPluginActSVG)
add_unittest(DD4hepDetectorElement DD4hepDetectorElementTests.cpp)
add_unittest(DD4hepCylinderLayerStructure DD4hepCylinderLayerStructureTests.cpp)
add_unittest(DD4hepDiscLayerStructure DD4hepDiscLayerStructureTests.cpp)

#set_property(TEST DD4hepDetectorElement PROPERTY ENVIRONMENT "PATH=path_list_append:/path/to/dlls")
#set_property(TEST DD4hepCylinderLayerStructure PROPERTY ENVIRONMENT "PATH=path_list_append:/path/to/dlls")
#set_property(TEST DD4hepDiscLayerStructure PROPERTY ENVIRONMENT "PATH=path_list_append:/path/to/dlls")
