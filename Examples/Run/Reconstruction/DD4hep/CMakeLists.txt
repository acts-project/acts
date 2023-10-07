add_executable(
  ActsExampleTruthTracksDD4hep
  DD4hepRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksDD4hep
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDD4hepDetectorWithOptions)

add_executable(
  ActsExampleCKFTracksDD4hep
  DD4hepRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksDD4hep
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDD4hepDetectorWithOptions)

add_executable(
  ActsExampleSeedingDD4hep
  DD4hepSeedingExample.cpp)
target_link_libraries(
  ActsExampleSeedingDD4hep
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDD4hepDetectorWithOptions)

install(
  TARGETS
    ActsExampleTruthTracksDD4hep
    ActsExampleCKFTracksDD4hep
    ActsExampleSeedingDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
