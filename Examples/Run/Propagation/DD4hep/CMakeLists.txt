add_executable(
  ActsExamplePropagationDD4hep
  DD4hepPropagationExample.cpp)
target_link_libraries(
  ActsExamplePropagationDD4hep
  PRIVATE ${_common_libraries} ActsExamplesDD4hepDetectorWithOptions)

install(
  TARGETS ActsExamplePropagationDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
