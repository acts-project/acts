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
    src/RootSpacePointPerformanceWriter.cpp
    src/RootSpacePointWriter.cpp
    src/RootTrackParameterPerformanceWriter.cpp
    src/RootTrackParameterWriter.cpp
    src/RootTrackStatesWriter.cpp
    src/RootTrackSummaryReader.cpp
    src/RootTrackSummaryWriter.cpp
    src/RootBFieldWriter.cpp
    src/RootAthenaNTupleReader.cpp
    src/RootVertexReader.cpp
    src/RootVertexWriter.cpp
    src/RootAthenaDumpReader.cpp
    src/RootAthenaDumpWriter.cpp
    src/RootNuclearInteractionParametersWriter.cpp
    src/detail/NuclearInteractionParametrisation.cpp
    src/RootTrackFinderPerformanceWriter.cpp
    src/RootTrackFinderNTupleWriter.cpp
    src/RootTrackFitterPerformanceWriter.cpp
    src/RootVertexNTupleWriter.cpp
    src/RootMuonSpacePointWriter.cpp
    src/RootMuonSpacePointReader.cpp
    src/RootMeasurementPerformanceWriter.cpp
    src/RootFileHasher.cpp
)
target_include_directories(
    ActsExamplesIoRoot
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesIoRoot
    PUBLIC
        Acts::PluginRoot
        Acts::ExamplesDigitization
        Acts::ExamplesFramework
        Acts::ExamplesPropagation
        Acts::ExamplesTruthTracking
        Acts::ExamplesMaterialMapping
        Threads::Threads
    PRIVATE ROOT::Core ROOT::Hist ROOT::Tree ROOT::TreePlayer
)

# RootFileHasher prefers Boost.Hash2 (SHA-256), which is available from Boost
# 1.86 onwards. On older Boost it falls back to ROOT's TMD5 (MD5).
if(Boost_VERSION VERSION_GREATER_EQUAL "1.86.0")
    target_compile_definitions(
        ActsExamplesIoRoot
        PRIVATE ACTS_ROOT_FILE_HASHER_USE_BOOST_HASH2
    )
    target_link_libraries(ActsExamplesIoRoot PRIVATE Boost::boost)
else()
    message(
        WARNING
        "Boost ${Boost_VERSION} < 1.86.0: RootFileHasher falls back to MD5 (ROOT TMD5) instead of SHA-256 (Boost.Hash2)."
    )
endif()

root_generate_dictionary(
  ActsExamplesIoRootDict MODULE ActsExamplesIoRoot LINKDEF LinkDef.hpp
)
set_target_properties(ActsExamplesIoRootDict PROPERTIES CXX_CLANG_TIDY "")
target_compile_options(ActsExamplesIoRootDict PRIVATE "-Wno-old-style-cast")
