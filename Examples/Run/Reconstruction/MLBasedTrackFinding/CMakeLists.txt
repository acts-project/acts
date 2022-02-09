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
    Boost::filesystem
)

install(
  TARGETS
  ActsExampleTrackFindingMLBased
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
