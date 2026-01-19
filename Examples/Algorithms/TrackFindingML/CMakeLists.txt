set(SOURCES
    src/AmbiguityResolutionMLAlgorithm.cpp
    src/SeedFilterMLAlgorithm.cpp
)

acts_add_library(ExamplesTrackFindingML ${SOURCES})

target_include_directories(
    ActsExamplesTrackFindingML
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTrackFindingML
    PUBLIC Acts::PluginOnnx Acts::ExamplesFramework
)

acts_compile_headers(ExamplesTrackFindingML GLOB "include/**/*.hpp")
