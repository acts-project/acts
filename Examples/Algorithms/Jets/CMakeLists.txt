acts_add_library(ExamplesJets SHARED src/TruthJetAlgorithm.cpp)
=======
add_library(ActsExamplesJets 
    SHARED 
    src/TruthJetAlgorithm.cpp
    src/TrackToTruthJetAlgorithm.cpp)
>>>>>>> ee1e83592 (Initialize TrackToTruthJetAlg)
target_include_directories(
    ActsExamplesJets
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesJets
    PUBLIC Acts::PluginFastJet Acts::ExamplesFramework
)

acts_compile_headers(ExamplesJets GLOB "include/**/*.hpp")
