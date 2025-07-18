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
    ExamplesTrackFindingExaTrkX
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
    ExamplesTrackFindingExaTrkX
    PUBLIC
        Acts::PluginExaTrkX
        Acts::ExamplesFramework
        Acts::ExamplesTrackFinding
)

install(
    TARGETS ActsExamplesTrackFindingExaTrkX
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
