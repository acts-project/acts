add_library(
	ActsExamplesRecTracksCommon SHARED
	RecTruthTracks.cpp
	RecCKFTracks.cpp)
target_include_directories(
	ActsExamplesRecTracksCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
	ActsExamplesRecTracksCommon
  PUBLIC
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFitting
    ActsExamplesTrackFinding
    ActsExamplesMagneticField
    ActsExamplesTruthTracking
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)    

install(
	TARGETS ActsExamplesRecTracksCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
