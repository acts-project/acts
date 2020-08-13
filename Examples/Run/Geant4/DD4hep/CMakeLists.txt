add_executable(
  ActsExampleDD4hepGeantinoRecording
  DD4hepGeantinoRecording.cpp)
target_link_libraries(
  ActsExampleDD4hepGeantinoRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4DD4hep
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsExampleDD4hepGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
