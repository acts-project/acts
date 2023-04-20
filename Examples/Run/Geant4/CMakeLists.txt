add_subdirectory(Common)

add_executable(
  ActsExampleMaterialRecordingGdml
  GdmlMaterialRecording.cpp)
target_link_libraries(
  ActsExampleMaterialRecordingGdml
  PRIVATE
    ActsExamplesCommon
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4Common
    ActsExamplesIoRoot)

add_executable(
  ActsExampleTestMockupBuilder
  TestMockupBuilder.cpp)
target_link_libraries(
  ActsExampleTestMockupBuilder
  PRIVATE
    ActsCore
    ActsExamplesCommon
    ActsExamplesFramework
    ActsExamplesGeant4
    ActsExamplesGeant4Common
    ActsExamplesIoRoot
    ActsPluginGeant4
    ActsExamplesDetectorGeant4
    ActsExamplesMuonSpectrometerMockupDetector
  )

install(
  TARGETS 
    ActsExampleMaterialRecordingGdml 
    ActsExampleTestMockupBuilder
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory_if(HepMC ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_HEPMC3)
