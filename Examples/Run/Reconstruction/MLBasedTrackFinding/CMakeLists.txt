add_executable(
  ActsExampleTrackFindingMLBased
  TrackFindingMLBasedExample.cpp)

target_link_libraries(
  ActsExampleTrackFindingMLBased
  PUBLIC
    ActsExamplesCommon
    ActsExamplesFramework
    ActsExamplesDigitization
    ActsExamplesDetectorGeneric
    ActsExamplesTrackFindingMLBased
    ActsPluginExaTrkX
    ActsExamplesIoPerformance
)

install(
  TARGETS
  ActsExampleTrackFindingMLBased
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})