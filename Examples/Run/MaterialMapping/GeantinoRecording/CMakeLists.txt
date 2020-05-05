add_executable(
  ActsExampleGeantinoRecording
  GeantinoRecordingExample.cpp)
target_include_directories(
  ActsExampleGeantinoRecording
  SYSTEM PRIVATE ${Geant4_INCLUDE_DIRS})
target_link_libraries(
  ActsExampleGeantinoRecording
  PRIVATE
    ${_common_libraries}
    ActsExamplesMaterialMapping ActsExamplesDD4hepG4
    ActsExamplesDetectorDD4hep ActsExamplesIoRoot
    ${Geant4_LIBRARIES})

install(
  TARGETS ActsExampleGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
