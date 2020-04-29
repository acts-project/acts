add_executable(
  ActsExampleGeometryDD4hep
  DD4hepGeometryExample.cpp)
target_include_directories(
  ActsExampleGeometryDD4hep
  PRIVATE ${DD4hep_INCLUDE_DIRS})
target_link_libraries(
  ActsExampleGeometryDD4hep
  PRIVATE ${_common_libraries} ActsExamplesDetectorDD4hep ${DD4hep_LIBRARIES})

install(
  TARGETS ActsExampleGeometryDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
