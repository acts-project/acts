acts_add_library(
    PluginRoot
    src/RootMaterialDecorator.cpp
    src/RootMaterialMapIo.cpp
    src/RootMaterialTrackIo.cpp
    src/RootMagneticFieldIo.cpp
    src/RootMeasurementIo.cpp
    src/RootSpacePointIo.cpp
    src/TGeoCylinderDiscSplitter.cpp
    src/TGeoDetectorElement.cpp
    src/TGeoLayerBuilder.cpp
    src/TGeoParser.cpp
    src/TGeoPrimitivesHelper.cpp
    src/TGeoSurfaceConverter.cpp
    src/HistogramConverter.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

# Internal system of units change, This behavior affected versions up to v6-25-01
# see: https://root.cern.ch/doc/v626/classTGeoMaterial.html
if(${ROOT_VERSION} VERSION_GREATER "6.25.01")
    target_sources(ActsPluginRoot PRIVATE src/TGeoMaterialConverter.cpp)
else()
    message(
        INFO
        " Skipping TGeoMaterialConverter due to ROOT version ${ROOT_VERSION}"
    )
endif()

target_include_directories(
    ActsPluginRoot
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
    ActsPluginRoot
    PUBLIC Acts::Core ROOT::Core ROOT::Tree ROOT::Geom
)

acts_compile_headers(PluginRoot GLOB include/**/*.hpp)
