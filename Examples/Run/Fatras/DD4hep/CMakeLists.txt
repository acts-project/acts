add_executable(
  ActsExampleFatrasDD4hep
  DD4hepFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasDD4hep
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorDD4hep)

install(
  TARGETS ActsExampleFatrasDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
