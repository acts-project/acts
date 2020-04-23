include(${Geant4_USE_FILE})

add_executable(
  ActsExampleGeantinoRecording
  GeantinoRecordingExample.cpp)
target_include_directories(
  ActsExampleGeantinoRecording
  PRIVATE ${DD4hep_INCLUDE_DIRS})
target_link_libraries(
  ActsExampleGeantinoRecording
  PRIVATE
    ${_common_libraries}
    ActsExamplesMaterialMapping ActsExamplesDD4hepG4
    ActsExamplesDetectorDD4hep ActsExamplesIoRoot
    ${DD4hep_LIBRARIES} ${Geant4_LIBRARIES})

install(
  TARGETS ActsExampleGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
