set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesPropagation
  ActsExamplesSeeding
  ActsExamplesCommon
  ActsExamplesMagneticField
  ActsExamplesIoObj
  ActsExamplesIoJson
  ActsExamplesIoRoot
  Boost::program_options)


add_executable(
  ActsExampleSeeding
  SeedingExample.cpp)

target_link_libraries(
  ActsExampleSeeding
  PRIVATE ${_common_libraries} ActsExamplesDetectorGeneric)


# install(
#   TARGETS
#     ActsExamplePropagationEmpty
#     ActsExamplePropagationGeneric
#     ActsExamplePropagationAligned
#     ActsExamplePropagationPayload
#     ActsExamplePropagationTGeo
#   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})


