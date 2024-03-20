add_executable(
  ActsExampleFatrasDD4hep
  DD4hepFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasDD4hep
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDD4hepDetectorWithOptions)

install(
  TARGETS ActsExampleFatrasDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
