add_library(
  ActsExamplesTrackFinding SHARED
  src/TrackFindingAlgorithm.cpp
  src/TrackFindingAlgorithmTrackFinderFunction.cpp)
target_include_directories(
  ActsExamplesTrackFinding
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTrackFinding
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField
    Boost::program_options)

install(
  TARGETS ActsExamplesTrackFinding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
