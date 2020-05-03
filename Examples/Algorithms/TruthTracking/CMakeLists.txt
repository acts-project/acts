add_library(
  ActsExamplesTruthTracking SHARED
  ACTFW/TruthTracking/ParticleSmearing.cpp
  ACTFW/TruthTracking/TrackSelector.cpp
  ACTFW/TruthTracking/TruthTrackFinder.cpp
  ACTFW/TruthTracking/TruthVerticesToTracks.cpp
  ACTFW/TruthTracking/TruthSeedSelector.cpp
  )
target_include_directories(
  ActsExamplesTruthTracking
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesTruthTracking
  PUBLIC ActsCore ActsExamplesFramework Boost::program_options)

install(
  TARGETS ActsExamplesTruthTracking
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
