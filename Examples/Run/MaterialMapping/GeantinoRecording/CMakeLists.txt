add_executable(
  ActsExampleGeantinoRecording
  GeantinoRecordingExample.cpp)
target_link_libraries(
  ActsExampleGeantinoRecording
  PRIVATE
    ${_common_libraries}
    ActsExamplesMaterialMapping ActsExamplesDD4hepG4
    ActsExamplesDetectorDD4hep ActsExamplesIoRoot)

install(
  TARGETS ActsExampleGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
