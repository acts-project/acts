acts_add_library(
    ExamplesIoRoot
    src/RootMeasurementWriter.cpp
    src/RootMaterialWriter.cpp
    src/RootMaterialTrackReader.cpp
    src/RootMaterialTrackWriter.cpp
    src/RootParticleWriter.cpp
    src/RootParticleReader.cpp
    src/RootPropagationStepsWriter.cpp
    src/RootPropagationSummaryWriter.cpp
    src/RootSeedWriter.cpp
    src/RootSimHitWriter.cpp
    src/RootSimHitReader.cpp
    src/RootSpacepointWriter.cpp
    src/RootTrackParameterWriter.cpp
    src/RootTrackStatesWriter.cpp
    src/RootTrackSummaryReader.cpp
    src/RootTrackSummaryWriter.cpp
    src/RootBFieldWriter.cpp
    src/RootAthenaNTupleReader.cpp
    src/RootVertexReader.cpp
    src/RootVertexWriter.cpp
    src/RootAthenaDumpReader.cpp
    src/RootNuclearInteractionParametersWriter.cpp
    src/detail/NuclearInteractionParametrisation.cpp
    src/RootTrackFinderPerformanceWriter.cpp
    src/RootTrackFinderNTupleWriter.cpp
    src/RootTrackFitterPerformanceWriter.cpp
    src/RootVertexNTupleWriter.cpp
    src/RootMuonSpacePointWriter.cpp
    src/RootMuonSpacePointReader.cpp
)
target_include_directories(
    ActsExamplesIoRoot
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoRoot
    PUBLIC
        Acts::PluginRoot
        Acts::ExamplesJets
        Acts::ExamplesDigitization
        Acts::ExamplesFramework
        Acts::ExamplesPropagation
        Acts::ExamplesTruthTracking
        Acts::ExamplesMaterialMapping
        Threads::Threads
    PRIVATE ROOT::Core ROOT::Hist ROOT::Tree ROOT::TreePlayer
)

root_generate_dictionary(
  ActsExamplesIoRootDict MODULE ActsExamplesIoRoot LINKDEF LinkDef.hpp
)
set_target_properties(ActsExamplesIoRootDict PROPERTIES CXX_CLANG_TIDY "")
target_compile_options(ActsExamplesIoRootDict PRIVATE "-Wno-old-style-cast")
