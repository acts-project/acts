add_library(
  ActsExamplesFramework SHARED
  src/Framework/BareAlgorithm.cpp
  src/Framework/BareService.cpp
  src/Framework/RandomNumbers.cpp
  src/Framework/Sequencer.cpp
  src/Utilities/Paths.cpp
  src/Utilities/Options.cpp
  src/Utilities/Helpers.cpp
  src/Validation/DuplicationPlotTool.cpp
  src/Validation/EffPlotTool.cpp
  src/Validation/FakeRatePlotTool.cpp
  src/Validation/ResPlotTool.cpp
  src/Validation/TrackClassification.cpp
  src/Validation/TrackSummaryPlotTool.cpp
)
target_include_directories(
  ActsExamplesFramework
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesFramework
  PUBLIC ActsCore ActsFatras Boost::boost ROOT::Core ROOT::Hist
  PRIVATE Boost::filesystem dfelibs)
target_compile_definitions(
  ActsExamplesFramework
  PRIVATE BOOST_FILESYSTEM_NO_DEPRECATED)
# set per-target c++17 requirement that will be propagated to linked targets
target_compile_features(
  ActsExamplesFramework
  PUBLIC cxx_std_17)

if(ACTS_USE_EXAMPLES_TBB)
  # newer DD4hep version require TBB and search interally for TBB in
  # config-only mode. to avoid mismatches we explicitely search using
  # config-only mode first to be sure that we find the same version.
  find_package(TBB ${_acts_tbb_version} CONFIG)
  if(NOT TBB_FOUND)
    # no version check possible when using the find module
    find_package(TBB ${_acts_tbb_version} MODULE)
  endif()
else()
  set(TBB_FOUND FALSE)
endif()
if(TBB_FOUND)
  target_link_libraries(
    ActsExamplesFramework
    PUBLIC TBB::tbb)
else()
  message(STATUS "disable TBB for Examples/Framework - only single-threaded running will be supported")
  target_compile_definitions(
    ActsExamplesFramework
    PUBLIC
    -DACTS_EXAMPLES_NO_TBB)
endif()

install(
  TARGETS ActsExamplesFramework
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
