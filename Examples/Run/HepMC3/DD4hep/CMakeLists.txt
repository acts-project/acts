add_executable(
  ActsExampleEventRecording
  EventRecordingExample.cpp)
target_link_libraries(
  ActsExampleEventRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesEventRecording
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4DD4hep
    ActsExamplesIoHepMC3
    Boost::program_options)
    
install(
  TARGETS ActsExampleEventRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
