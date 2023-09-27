set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesPropagation
  ActsExamplesCommon
  ActsExamplesMagneticField
  ActsExamplesIoObj
  ActsExamplesIoJson
  ActsExamplesIoRoot)

# Empty detector
add_executable(
  ActsExamplePropagationEmpty
  EmptyDetectorPropagationExample.cpp)
target_link_libraries(
  ActsExamplePropagationEmpty
  PRIVATE ${_common_libraries})

# Generic detector
add_executable(
  ActsExamplePropagationGeneric
  GenericPropagationExample.cpp)
target_link_libraries(
  ActsExamplePropagationGeneric
  PRIVATE ${_common_libraries})

# Generic detector with alignment
add_executable(
  ActsExamplePropagationAligned
  AlignedPropagationExample.cpp)
target_link_libraries(
  ActsExamplePropagationAligned
  PRIVATE ${_common_libraries})


# TGeo detector
add_executable(
  ActsExamplePropagationTGeo
  TGeoPropagationExample.cpp)
target_link_libraries(
  ActsExamplePropagationTGeo
  PRIVATE ${_common_libraries})

install(
  TARGETS
    ActsExamplePropagationEmpty
    ActsExamplePropagationGeneric
    ActsExamplePropagationAligned
    ActsExamplePropagationTGeo
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
