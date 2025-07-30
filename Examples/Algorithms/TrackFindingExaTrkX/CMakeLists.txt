acts_add_library(
    ExamplesTrackFindingExaTrkX
    SHARED
    src/TrackFindingAlgorithmExaTrkX.cpp
    src/PrototracksToParameters.cpp
    src/TrackFindingFromPrototrackAlgorithm.cpp
    src/TruthGraphBuilder.cpp
    src/createFeatures.cpp
)

target_include_directories(
    ActsExamplesTrackFindingExaTrkX
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTrackFindingExaTrkX
    PUBLIC
        Acts::PluginExaTrkX
        Acts::ExamplesFramework
        Acts::ExamplesTrackFinding
)
