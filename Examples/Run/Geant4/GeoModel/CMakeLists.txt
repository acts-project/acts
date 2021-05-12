add_executable(
  ActsExampleGeantinoRecordingGdml
  GdmlGeantinoRecording.cpp)
target_link_libraries(
  ActsExampleGeantinoRecordingGdml
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4GeoModel
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS 
    ActsExampleGeantinoRecordingGdml 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})