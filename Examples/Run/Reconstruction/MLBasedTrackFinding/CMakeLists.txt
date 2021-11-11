add_executable(
  ActsExampleTrackFindingMLBased
  TrackFindingMLBasedExample.cpp)

target_link_libraries(
  ActsExampleTrackFindingMLBased
  PRIVATE
    ActsExamplesCommon
    ActsExamplesDigitization
    ActsExamplesDetectorGeneric
    ActsExamplesTrackFindingMLBased
)

install(
  TARGETS
  ActsExampleTrackFindingMLBased
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})