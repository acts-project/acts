acts_add_library(
    ExamplesTrackFindingGnn
    src/TrackFindingAlgorithmGnn.cpp
    src/ProtoTracksToParameters.cpp
    src/TrackFindingFromProtoTrackAlgorithm.cpp
    src/TruthGraphBuilder.cpp
    src/createFeatures.cpp
)

target_include_directories(
    ActsExamplesTrackFindingGnn
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ActsExamplesTrackFindingGnn
    PUBLIC Acts::PluginGnn Acts::ExamplesFramework Acts::ExamplesTrackFinding
)
