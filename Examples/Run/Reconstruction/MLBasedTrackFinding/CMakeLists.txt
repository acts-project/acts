add_executable(
  ActsExampleTrackFindingMLBased
  TrackFindingMLBasedExample.cpp)

target_link_libraries(
  ActsExampleTrackFindingMLBased
  PRIVATE
    ActsExamplesCommon
    ActsExamplesFramework
    ActsExamplesDigitization
    ActsExamplesDetectorGeneric
    ActsExamplesIoPerformance
    ActsExamplesTrackFindingMLBased
    ActsPluginExaTrkX
)

install(
  TARGETS
  ActsExampleTrackFindingMLBased
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})