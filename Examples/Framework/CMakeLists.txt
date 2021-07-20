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
  PUBLIC ActsCore ActsFatras Boost::boost ROOT::Core ROOT::Hist TBB::tbb
  PRIVATE Boost::filesystem dfelibs)
target_compile_definitions(
  ActsExamplesFramework
  PRIVATE BOOST_FILESYSTEM_NO_DEPRECATED)
# set per-target c++17 requirement that will be propagated to linked targets
target_compile_features(
  ActsExamplesFramework
  PUBLIC cxx_std_17)

install(
  TARGETS ActsExamplesFramework
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# if(ACTS_BUILD_PLUGIN_PYTHON)
# pybind11_add_module()
# endif()