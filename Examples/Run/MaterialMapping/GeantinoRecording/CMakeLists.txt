add_executable(
  ActsExampleGeantinoRecording
  GeantinoRecordingExample.cpp)
target_link_libraries(
  ActsExampleGeantinoRecording
  PRIVATE
    ${_common_libraries}
    ActsExamplesDetectorDD4hep
    ActsExamplesGeant4 ActsExamplesIoRoot ActsExamplesMaterialMapping)

install(
  TARGETS ActsExampleGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
