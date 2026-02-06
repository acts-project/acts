include(ActsTargetLinkLibrariesSystem)

acts_add_library(
    ExamplesFramework
    src/EventData/MuonSpacePoint.cpp
    src/EventData/MuonSpacePointCalibrator.cpp
    src/EventData/Measurement.cpp
    src/EventData/MeasurementCalibration.cpp
    src/EventData/SimParticle.cpp
    src/Framework/IAlgorithm.cpp
    src/Framework/SequenceElement.cpp
    src/Framework/WhiteBoard.cpp
    src/Framework/RandomNumbers.cpp
    src/Framework/Sequencer.cpp
    src/Framework/DataHandle.cpp
    src/Framework/BufferedReader.cpp
    src/Utilities/EventDataTransforms.cpp
    src/Utilities/Paths.cpp
    src/Utilities/Options.cpp
    src/Utilities/ParametricParticleGenerator.cpp
    src/Validation/DuplicationPlotTool.cpp
    src/Validation/EffPlotTool.cpp
    src/Validation/FakePlotTool.cpp
    src/Validation/ResPlotTool.cpp
    src/Validation/TrackClassification.cpp
    src/Validation/TrackQualityPlotTool.cpp
    src/Validation/TrackSummaryPlotTool.cpp
    ACTS_INCLUDE_FOLDER include/ActsExamples
)

target_include_directories(
    ActsExamplesFramework
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(
    ActsExamplesFramework
    PUBLIC Acts::Fatras Acts::PluginFpeMonitoring Boost::boost HepMC3::HepMC3
)
target_compile_definitions(
    ActsExamplesFramework
    PRIVATE BOOST_FILESYSTEM_NO_DEPRECATED
)

# Plugin-dependent framework components
# Note that due to conditional compilation the headers of ROOT/Onnx components
# get installed regardless of the configuration which could lead to link errors
# in case one includes these files without building the components.
if(ACTS_BUILD_EXAMPLES_ROOT)
    target_sources(
        ActsExamplesFramework
        PRIVATE src/Root/ScalingCalibrator.cpp src/Root/MuonVisualization.cpp
    )
    target_link_libraries(
        ActsExamplesFramework
        PUBLIC ROOT::Core ROOT::Hist ROOT::Graf ROOT::Gpad
    )
    set_source_files_properties(
        src/Framework/Sequencer.cpp
        TARGET
        ActsExamplesFramework
        PROPERTIES COMPILE_DEFINITIONS ACTS_BUILD_EXAMPLES_ROOT
    )
endif()

if(ACTS_BUILD_PLUGIN_ONNX)
    target_sources(ActsExamplesFramework PRIVATE src/Onnx/NeuralCalibrator.cpp)
    target_link_libraries(ActsExamplesFramework PUBLIC Acts::PluginOnnx)
endif()

# newer DD4hep version require TBB and search internally for TBB in
# config-only mode. to avoid mismatches we explicitly search using
# config-only mode first to be sure that we find the same version.
find_package(TBB ${_acts_tbb_version} CONFIG)

if(NOT TBB_FOUND)
    # no version check possible when using the find module
    find_package(TBB ${_acts_tbb_version} MODULE REQUIRED)
endif()

target_link_libraries(ActsExamplesFramework PUBLIC TBB::tbb)
