acts_add_library(
    ExamplesJets
    SHARED
    src/TruthJetAlgorithm.cpp
    src/TrackToTruthJetAlgorithm.cpp
)
target_include_directories(
    ActsExamplesJets
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesJets
    PUBLIC ActsCore ActsPluginFastJet ActsExamplesFramework Acts::ExamplesIoHepMC3 HepMC3::HepMC3
)

acts_compile_headers(ExamplesJets GLOB "include/**/*.hpp")
