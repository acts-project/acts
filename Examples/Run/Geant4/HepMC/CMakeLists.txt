add_executable(
  ActsExampleEventRecording
  EventRecordingExample.cpp)
target_link_libraries(
  ActsExampleEventRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDD4hepDetectorWithOptions
    ActsExamplesEventRecording
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesIoHepMC3
    ActsExamplesHepMC3Options
    ActsExamplesGeant4Common)
    
install(
  TARGETS ActsExampleEventRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
