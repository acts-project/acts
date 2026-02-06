acts_add_library(
    ExamplesAmbiguityResolution
    src/GreedyAmbiguityResolutionAlgorithm.cpp
    src/ScoreBasedAmbiguityResolutionAlgorithm.cpp
)

target_include_directories(
    ActsExamplesAmbiguityResolution
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesAmbiguityResolution
    PUBLIC Acts::ExamplesFramework Acts::PluginJson
)
