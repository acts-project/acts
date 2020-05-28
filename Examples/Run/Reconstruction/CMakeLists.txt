add_executable(
  ActsRecTruthTracks
  RecTruthTracks.cpp)

target_link_libraries(
  ActsRecTruthTracks
  PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesFitting
    ActsExamplesDetectorGeneric
    ActsExamplesMagneticField
    ActsExamplesTruthTracking
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS ActsRecTruthTracks
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
    ActsRecCKFTracks
    RecCKFTracks.cpp)
target_link_libraries(
    ActsRecCKFTracks
    PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesTrackFinding
    ActsExamplesDetectorGeneric
    ActsExamplesMagneticField
    ActsExamplesIoCsv
    ActsExamplesIoPerformance)

install(
  TARGETS
  ActsRecCKFTracks
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

