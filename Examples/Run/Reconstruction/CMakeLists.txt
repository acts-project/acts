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
