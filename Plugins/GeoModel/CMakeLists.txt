acts_add_library(
    PluginGeoModel
    src/GeoModelTree.cpp
    src/GeoModelConversionError.cpp
    src/GeoModelToDetectorVolume.cpp
    src/GeoModelReader.cpp
    src/GeoModelDetectorElement.cpp
    src/GeoModelDetectorObjectFactory.cpp
    src/GeoModelMaterialConverter.cpp
    src/detail/GeoBoxConverter.cpp
    src/detail/GeoTrdConverter.cpp
    src/detail/GeoTubeConverter.cpp
    src/detail/GeoPolygonConverter.cpp
    src/detail/GeoSubtractionConverter.cpp
    src/detail/GeoShiftConverter.cpp
    src/detail/GeoIntersectionAnnulusConverter.cpp
    src/detail/GeoModelExtentHelper.cpp
    src/detail/GeoUnionDoubleTrdConverter.cpp
    src/GeoModelDetectorElementITk.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginGeoModel
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginGeoModel
    PUBLIC
        Acts::Core
        GeoModelCore::GeoModelKernel
        GeoModelIO::GeoModelDBManager
        GeoModelIO::GeoModelRead
)

if(ACTS_BUILD_PLUGIN_JSON)
    target_link_libraries(ActsPluginGeoModel PRIVATE ActsPluginJson)
    target_sources(
        ActsPluginGeoModel
        PRIVATE src/GeoModelJsonMaterialManager.cpp
    )
    target_compile_definitions(
        ActsPluginGeoModel
        PUBLIC ACTS_BUILD_PLUGIN_JSON=1
    )
endif()

acts_compile_headers(PluginGeoModel GLOB include/**/*.hpp)
