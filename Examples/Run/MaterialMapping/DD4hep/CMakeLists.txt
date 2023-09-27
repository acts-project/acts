add_executable(
  ActsExampleMaterialValidationDD4hep
  DD4hepMaterialValidation.cpp)
target_link_libraries(
  ActsExampleMaterialValidationDD4hep
  PRIVATE
    ${_common_libraries}
    ActsExamplesPropagation ActsExamplesDD4hepDetectorWithOptions)

add_executable(
  ActsExampleMaterialMappingDD4hep
  DD4hepMaterialMapping.cpp)
target_link_libraries(
  ActsExampleMaterialMappingDD4hep
  PRIVATE
    ${_common_libraries}
    ActsExamplesMaterialMapping ActsExamplesDD4hepDetectorWithOptions)

install(
  TARGETS
    ActsExampleMaterialValidationDD4hep
    ActsExampleMaterialMappingDD4hep
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
