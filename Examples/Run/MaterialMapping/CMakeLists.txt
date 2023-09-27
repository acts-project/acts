set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesIoRoot
  ActsExamplesIoJson
  ActsExamplesMagneticField)

add_executable(
  ActsExampleMaterialValidationGeneric
  GenericMaterialValidation.cpp)
target_link_libraries(
  ActsExampleMaterialValidationGeneric
  PRIVATE ${_common_libraries} ActsExamplesPropagation ActsExamplesDetectorGeneric)

add_executable(
  ActsExampleMaterialMappingGeneric
  GenericMaterialMapping.cpp)
target_link_libraries(
  ActsExampleMaterialMappingGeneric
  PRIVATE ${_common_libraries} ActsExamplesMaterialMapping ActsExamplesDetectorGeneric)

install(
  TARGETS
    ActsExampleMaterialValidationGeneric
    ActsExampleMaterialMappingGeneric
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
add_subdirectory_if(TGeo ACTS_BUILD_PLUGIN_TGEO)
