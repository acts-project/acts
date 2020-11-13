add_library(
  ActsExamplesRecTracksCommon SHARED
  RecTruthTracks.cpp
  RecCKFTracks.cpp)
target_include_directories(
  ActsExamplesRecTracksCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
if(ACTS_BUILD_PLUGIN_ONNX)
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
      ActsExamplesIoPerformance
      ActsPluginOnnx)
else()
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
endif()

install(
  TARGETS ActsExamplesRecTracksCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
