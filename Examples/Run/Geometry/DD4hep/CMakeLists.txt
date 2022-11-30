add_executable(
  ActsExampleGeometryDD4hep
  DD4hepGeometryExample.cpp)
target_link_libraries(
  ActsExampleGeometryDD4hep
  PRIVATE ${_common_libraries} ActsExamplesDD4hepDetectorWithOptions)

install(
  TARGETS ActsExampleGeometryDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
