add_executable(
  ActsExampleTruthTracks
  RecTruthTracks.cpp)

target_link_libraries(
  ActsExampleTruthTracks
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
  TARGETS ActsExampleTruthTracks
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(
    ActsExampleCKFTracks
    RecCKFTracks.cpp)
target_link_libraries(
    ActsExampleCKFTracks
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
  ActsExampleCKFTracks
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

