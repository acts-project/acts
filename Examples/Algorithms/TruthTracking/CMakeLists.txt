acts_add_library(
    ExamplesTruthTracking
    ActsExamples/TruthTracking/ParticleSelector.cpp
    ActsExamples/TruthTracking/ParticleTrackParamExtractor.cpp
    ActsExamples/TruthTracking/TrackParameterSmearing.cpp
    ActsExamples/TruthTracking/TrackParameterSelector.cpp
    ActsExamples/TruthTracking/TrackModifier.cpp
    ActsExamples/TruthTracking/TrackTruthMatcher.cpp
    ActsExamples/TruthTracking/TruthTrackFinder.cpp
    ActsExamples/TruthTracking/TruthVertexFinder.cpp
    ActsExamples/TruthTracking/TruthSeedingAlgorithm.cpp
    ActsExamples/TruthTracking/HitSelector.cpp
)
target_include_directories(
    ActsExamplesTruthTracking
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(ActsExamplesTruthTracking PUBLIC Acts::ExamplesFramework)
