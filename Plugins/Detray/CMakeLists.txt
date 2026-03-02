acts_add_library(
    PluginDetray
    src/DetrayConversionUtils.cpp
    #src/DetrayConverter.cpp
    #src/DetrayGeometryConverter.cpp
    #src/DetrayMaterialConverter.cpp
    #src/DetraySurfaceGridsConverter.cpp
    src/DetrayPayloadConverter.cpp
    src/DetrayMaterial.cpp
    src/DetrayNavigation.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

add_dependencies(ActsPluginDetray detray::core covfie::core vecmem::core)

target_include_directories(
    ActsPluginDetray
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    ActsPluginDetray
    PUBLIC
        Acts::PluginJson
        detray::core
        detray::core_array
        detray::io
        detray::detectors
        detray::test_utils
        detray::svgtools
        vecmem::core
)

acts_compile_headers(PluginDetray GLOB include/**/*.hpp)
