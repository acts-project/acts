add_executable(
  ActsSimGeantinoRecording
  RunSimGeantinoRecording.cpp)
target_link_libraries(
  ActsSimGeantinoRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4DD4hep
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsSimGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
