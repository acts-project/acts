acts_add_library(ExamplesDetray INTERFACE)

target_include_directories(
    ActsExamplesDetray
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesDetray
    INTERFACE
        Acts::ExamplesFramework
        Acts::ExamplesPropagation
        Acts::PluginDetray
        Acts::PluginCovfie
)

acts_compile_headers(ExamplesDetray GLOB "include/**/*.hpp")
