acts_add_library(
    PluginActSVG
    SHARED
    src/DetectorVolumeSvgConverter.cpp
    src/DetectorSvgConverter.cpp
    src/EventDataSvgConverter.cpp
    src/LayerSvgConverter.cpp
    src/PortalSvgConverter.cpp
    src/SurfaceArraySvgConverter.cpp
    src/SurfaceSvgConverter.cpp
    src/TrackingGeometrySvgConverter.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginActSVG
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(
    ActsPluginActSVG
    PUBLIC Acts::Core actsvg::core actsvg::meta
)
