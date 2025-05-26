add_library(ActsExamplesJets SHARED src/TruthJetAlgorithm.cpp)
target_include_directories(
    ActsExamplesJets
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesJets
    PUBLIC ActsCore ActsPluginFastJet ActsExamplesFramework
)

acts_compile_headers(ActsExamplesJets GLOB "include/**/*.hpp")

install(TARGETS ActsExamplesJets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
