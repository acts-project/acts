add_executable(
  ActsExampleTruthTracksDD4hep
  DD4hepRecTruthTracks.cpp)

target_link_libraries(
  ActsExampleTruthTracksDD4hep
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorDD4hep)

add_executable(
  ActsExampleCKFTracksDD4hep
  DD4hepRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksDD4hep
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorDD4hep)

install(
  TARGETS
    ActsExampleTruthTracksDD4hep
    ActsExampleCKFTracksDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
