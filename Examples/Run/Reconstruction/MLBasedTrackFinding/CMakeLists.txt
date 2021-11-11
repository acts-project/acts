add_executable(
  ActsExampleTrackFindingMLBased
  TrackFindingMLBasedExample.cpp)
target_link_libraries(
  ActsExampleTrackFindingMLBased
  PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesIoCsv
    ActsExamplesIoJson
    ActsExamplesIoRoot
    ActsExamplesDetectorGeneric
    ActsExamplesDetectorTGeo
    ActsExamplesDigitization
    ActsExamplesMagneticField
    ActsExamplesTrackFinding # spacepoint maker should not live here!
    ActsExamplesTrackFindingMLBased
    ActsExamplesTrackFitting
  )


install(
  TARGETS
  ActsExampleTrackFindingMLBased
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})