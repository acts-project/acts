acts_add_library(ExamplesJets src/TruthJetAlgorithm.cpp)
target_include_directories(
    ActsExamplesJets
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesJets
    PUBLIC Acts::Fatras Acts::PluginFastJet Acts::ExamplesFramework
)

acts_compile_headers(ExamplesJets GLOB "include/**/*.hpp")
