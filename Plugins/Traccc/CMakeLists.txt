acts_add_library(
    PluginTraccc
    src/TrackConversion.cpp
    ACTS_INCLUDE_FOLDER include/ActsPlugins
)

target_include_directories(
    ActsPluginTraccc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(
    ActsPluginTraccc
    PUBLIC Acts::Core traccc::core detray::core vecmem::core
)

acts_compile_headers(PluginTraccc GLOB include/**/*.hpp)
