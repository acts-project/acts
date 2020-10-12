add_executable(
  ActsExampleHepMC3
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3
  PRIVATE
    ActsCore ActsExamplesCommon
    ActsExamplesFramework ActsExamplesIoHepMC3
    Boost::program_options)

add_executable(
  ActsExampleEventRecording
  EventRecordingExample.cpp)
target_link_libraries(
  ActsExampleEventRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4DD4hep
    ActsExamplesIoHepMC3
    Boost::program_options)
    
install(
  TARGETS ActsExampleHepMC3 ActsExampleEventRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
