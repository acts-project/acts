add_executable(
  ActsExamplePropagationDD4hep
  DD4hepPropagationExample.cpp)
target_link_libraries(
  ActsExamplePropagationDD4hep
  PRIVATE ${_common_libraries} ActsExamplesDetectorDD4hep)

install(
  TARGETS ActsExamplePropagationDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
