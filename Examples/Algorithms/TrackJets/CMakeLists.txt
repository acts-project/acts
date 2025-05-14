add_library(
    ActsExamplesTrackJets
    SHARED
    src/TruthJetAlgorithm.cpp
)
target_include_directories(
    ActsExamplesTrackJets
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesTrackJets
    PUBLIC ActsCore ActsPluginFastJet ActsExamplesFramework
)

acts_compile_headers(ActsExamplesTrackJets GLOB "include/**/*.hpp")

install(
    TARGETS ActsExamplesTrackJets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
