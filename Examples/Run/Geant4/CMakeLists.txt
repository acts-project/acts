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
    ActsExamplesIoRoot
    Boost::program_options)

add_executable(
  ActsNuclearInteractionRecording
  NuclearInteractionRecording.cpp)
target_link_libraries(
  ActsNuclearInteractionRecording
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDetectorDD4hep
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4DD4hep
    ActsExamplesIoHepMC3
    Boost::program_options)
    
install(
  TARGETS 
    ActsExampleGeantinoRecordingGdml 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory_if(HepMC ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_HEPMC3)
