acts_add_library(
    PluginDetray
    src/DetrayConversionUtils.cpp
    src/DetrayPayloadConverter.cpp
    src/DetrayMaterial.cpp
    src/DetrayNavigation.cpp
    src/DetrayDetectorOdd.cpp
    src/DetrayDetectorDefault.cpp
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
        detray::svgtools
        detray::test_utils
        vecmem::core
)

acts_compile_headers(PluginDetray GLOB include/**/*.hpp)
