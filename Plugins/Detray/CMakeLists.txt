acts_add_library(
    PluginDetray
    src/DetrayConversionUtils.cpp
    src/DetrayPayloadConverter.cpp
    src/DetrayMaterial.cpp
    src/DetrayNavigation.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

# The heavy detray detector operations (detector assembly, consistency check,
# JSON read/write) are explicitly instantiated for the closed set of shipped
# metadata types inside detray itself (see detray::detector_io_array); this
# plugin only links the pre-built symbols and sees the matching `extern
# template` declarations via detray's detector_io_array.hpp.

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
        detray::detector_io_array
        detray::svgtools
        detray::test_utils
        vecmem::core
)

acts_compile_headers(PluginDetray GLOB include/**/*.hpp)
