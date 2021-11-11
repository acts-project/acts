add_executable(
  ActsExampleTrackFindingMLBased
  TrackFindingMLBasedExample.cpp)

target_link_libraries(
  ActsExampleTrackFindingMLBased
  PRIVATE
    # ActsExamplesCommon
    # ActsExamplesFramework
    # ActsExamplesDigitization
    # ActsExamplesDetectorGeneric
    # ActsExamplesTrackFindingMLBased
    ActsExamplesTrackFinding
    # frnn cugraphtestutil ${TORCH_LIBRARIES} ${OnnxRuntime_LIBRARY}
)

install(
  TARGETS
  ActsExampleTrackFindingMLBased
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})