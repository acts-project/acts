acts_add_library(
    PluginDD4hep
    src/GeometryModuleHelper.cpp
    src/GeometryModuleLoader.cpp
    src/ConvertDD4hepDetector.cpp
    src/DD4hepBinningHelpers.cpp
    src/DD4hepMaterialHelpers.cpp
    src/DD4hepDetectorElement.cpp
    src/DD4hepDetectorSurfaceFactory.cpp
    src/DD4hepLayerBuilder.cpp
    src/DD4hepVolumeBuilder.cpp
    src/DD4hepFieldAdapter.cpp
    src/BlueprintBuilder.cpp
    src/OpenDataDetectorBuilder.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

set(_acts_dd4hep_geometry_module_abi_tag
    "${ACTS_GEOMETRY_MODULE_ABI_TAG}|dd4hep-${DD4hep_VERSION}"
)

target_include_directories(
    ActsPluginDD4hep
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(ActsPluginDD4hep PUBLIC Acts::PluginRoot)
target_compile_definitions(
    ActsPluginDD4hep
    PRIVATE
        ACTS_DD4HEP_GEOMETRY_MODULE_ABI_TAG="${_acts_dd4hep_geometry_module_abi_tag}"
)

acts_compile_headers(PluginDD4hep GLOB include/**/*.hpp)

if(${DD4hep_VERSION} VERSION_LESS 1.11)
    target_include_directories(ActsPluginDD4hep PUBLIC ${DD4hep_INCLUDE_DIRS})
    target_link_libraries(ActsPluginDD4hep PUBLIC ${DD4hep_DDCORE_LIBRARY})
else()
    target_link_libraries(ActsPluginDD4hep PUBLIC DD4hep::DDCore DD4hep::DDRec)
endif()
