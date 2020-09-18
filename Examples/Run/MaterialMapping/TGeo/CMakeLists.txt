add_executable(
  ActsExampleMaterialValidationTGeo
  TGeoMaterialValidation.cpp)
target_link_libraries(
  ActsExampleMaterialValidationTGeo
  PRIVATE
    ${_common_libraries}
    ActsExamplesPropagation ActsExamplesDetectorTGeo)

add_executable(
  ActsExampleMaterialMappingTGeo
  TGeoMaterialMapping.cpp)
target_link_libraries(
  ActsExampleMaterialMappingTGeo
  PRIVATE
    ${_common_libraries}
    ActsExamplesMaterialMapping ActsExamplesDetectorTGeo)

install(
  TARGETS
    ActsExampleMaterialValidationTGeo
    ActsExampleMaterialMappingTGeo
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
