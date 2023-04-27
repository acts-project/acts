add_library(
  ActsExamplesTrackFittingChi2 SHARED
  src/TrackFittingChi2Algorithm.cpp
  src/TrackFittingChi2AlgorithmFunction.cpp)
target_include_directories(
  ActsExamplesTrackFittingChi2
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTrackFittingChi2
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField)

install(
  TARGETS ActsExamplesTrackFittingChi2
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
