add_library(
  ActsExamplesTrackFinding SHARED
  src/AmbiguityResolutionAlgorithm.cpp
  src/SeedingAlgorithm.cpp
  src/SeedingOrthogonalAlgorithm.cpp
  src/SpacePointMaker.cpp
  src/TrackFindingAlgorithm.cpp
  src/TrackFindingAlgorithmFunction.cpp
  src/HoughTransformSeeder.cpp
  src/TrackParamsEstimationAlgorithm.cpp
  src/TrajectoriesToPrototracks.cpp)
target_include_directories(
  ActsExamplesTrackFinding
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTrackFinding
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ActsExamplesIoJson
    ActsExamplesMagneticField)

install(
  TARGETS ActsExamplesTrackFinding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
