add_executable(
  ActsExampleExaTrkXGeneric
  TrackFindingExaTrkXGeneric.cpp)

target_link_libraries(
  ActsExampleExaTrkXGeneric
  PRIVATE
    ActsExamplesCommon
    ActsExamplesFramework
    ActsExamplesDigitization
    ActsExamplesDetectorGeneric
    ActsExamplesIoPerformance
    ActsExamplesTrackFindingExaTrkX
    ActsPluginExaTrkX
    Boost::filesystem
)

install(
  TARGETS
  ActsExampleExaTrkXGeneric
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
