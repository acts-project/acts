add_executable(
  ActsExampleGdmlGeantinoRecording
  GdmlGeantinoRecording.cpp)
target_link_libraries(
  ActsExampleGdmlGeantinoRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesIoRoot
    Boost::program_options)

install(
  TARGETS ActsExampleGdmlGeantinoRecording
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
