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

# traccc is linked privately, so that the exported plugin target does not
# reference traccc::core: the Acts package configuration does not look up
# traccc (see ActsConfig.cmake.in). Users of the plugin headers have to find
# and link traccc themselves.
target_link_libraries(
    ActsPluginTraccc
    PUBLIC Acts::Core detray::core vecmem::core
    PRIVATE traccc::core
)

acts_compile_headers(PluginTraccc GLOB include/**/*.hpp)

# The header check needs traccc, which is not part of the plugin interface
if(TARGET ActsPluginTraccc_HEADERS)
    target_link_libraries(ActsPluginTraccc_HEADERS PRIVATE traccc::core)
endif()
