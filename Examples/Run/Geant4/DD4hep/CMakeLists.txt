add_executable(
  ActsExampleGeantinoRecordingDD4hep
  DD4hepGeantinoRecording.cpp)
target_link_libraries(
  ActsExampleGeantinoRecordingDD4hep
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4DD4hep
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsExampleGeantinoRecordingDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
