add_executable(
  ActsExamplePropagationDD4hep
  DD4hepPropagationExample.cpp)
target_include_directories(
  ActsExamplePropagationDD4hep
  PRIVATE ${DD4hep_INCLUDE_DIRS})
target_link_libraries(
  ActsExamplePropagationDD4hep
  PRIVATE ${_common_libraries} ActsExamplesDetectorDD4hep ${DD4hep_LIBRARIES})

install(
  TARGETS ActsExamplePropagationDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
