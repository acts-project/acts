add_executable(
  ActsSimFatrasDD4hep
  DD4hepFatrasExample.cpp)
target_include_directories(
  ActsSimFatrasDD4hep
  PRIVATE ${DD4hep_INCLUDE_DIRS})
target_link_libraries(
  ActsSimFatrasDD4hep
  PRIVATE
    ActsExamplesFatrasCommon ActsExamplesDetectorDD4hep
    ${DD4hep_LIBRARIES})

install(
  TARGETS ActsSimFatrasDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
