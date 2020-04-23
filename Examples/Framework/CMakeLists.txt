add_library(
  ActsExamplesFramework SHARED
  src/Framework/BareAlgorithm.cpp
  src/Framework/BareService.cpp
  src/Framework/RandomNumbers.cpp
  src/Framework/Sequencer.cpp
  src/Utilities/Paths.cpp
  src/Utilities/Options.cpp
  src/Utilities/Helpers.cpp
  src/Validation/EffPlotTool.cpp
  src/Validation/FakeRatePlotTool.cpp
  src/Validation/TrackSummaryPlotTool.cpp
  src/Validation/ProtoTrackClassification.cpp
  src/Validation/ResPlotTool.cpp)
target_include_directories(
  ActsExamplesFramework
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(
  ActsExamplesFramework
  PUBLIC ActsCore ActsFatras Boost::boost ROOT::Core ROOT::Hist
  PRIVATE ${TBB_LIBRARIES} Boost::filesystem dfelibs)
target_compile_definitions(
  ActsExamplesFramework
  # framework event data is based on a specific identifier
  PUBLIC ACTS_CORE_IDENTIFIER_PLUGIN="${CMAKE_CURRENT_SOURCE_DIR}/include/ACTFW/EventData/SimIdentifier.hpp"
  PRIVATE BOOST_FILESYSTEM_NO_DEPRECATED)
# set per-target c++17 requirement that will be propagated to linked targets
target_compile_features(
  ActsExamplesFramework
  PUBLIC cxx_std_17)

install(
  TARGETS ActsExamplesFramework
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
