add_executable(
  ActsSimFatrasDD4hep
  DD4hepFatrasExample.cpp)
target_link_libraries(
  ActsSimFatrasDD4hep
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorDD4hep)

install(
  TARGETS ActsSimFatrasDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
