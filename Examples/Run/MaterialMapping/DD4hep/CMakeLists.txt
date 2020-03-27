add_executable(
  ActsExampleMaterialValidationDD4hep
  DD4hepMaterialValidation.cpp)
target_link_libraries(
  ActsExampleMaterialValidationDD4hep
  PRIVATE
    ${_common_libraries}
    ActsExamplesPropagation ActsExamplesDetectorDD4hep
    ${DD4hep_LIBRARIES})

add_executable(
  ActsExampleMaterialMappingDD4hep
  DD4hepMaterialMapping.cpp)
target_link_libraries(
  ActsExampleMaterialMappingDD4hep
  PRIVATE
    ${_common_libraries}
    ActsExamplesMaterialMapping ActsExamplesDetectorDD4hep
    ${DD4hep_LIBRARIES})

install(
  TARGETS
    ActsExampleMaterialValidationDD4hep
    ActsExampleMaterialMappingDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
