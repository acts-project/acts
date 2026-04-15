acts_add_library(ExamplesTraccc INTERFACE)

target_include_directories(
    ActsExamplesTraccc
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTraccc
    INTERFACE
        Acts::ExamplesFramework
        Acts::ExamplesPropagation
        Acts::PluginDetray
        Acts::PluginCovfie
        traccc::io
        traccc::core
        $<$<BOOL:${ACTS_BUILD_PLUGIN_TRACCC}>:traccc::cuda>
        $<$<BOOL:${ACTS_BUILD_PLUGIN_TRACCC}>:vecmem::cuda>
)

acts_compile_headers(ExamplesTraccc GLOB "include/**/*.hpp")
