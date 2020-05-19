# Common library that implements the core functionality
add_subdirectory(Common)

# Generic detector
add_executable(
  ActsSimFatrasGeneric
  GenericFatrasExample.cpp)
target_link_libraries(
  ActsSimFatrasGeneric
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorGeneric)

# Generic detector with IOV based alignment
add_executable(
  ActsSimFatrasAligned
  AlignedFatrasExample.cpp)
target_link_libraries(
  ActsSimFatrasAligned
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorContextual)

# Generic detector with Payload based alignment
add_executable(
  ActsSimFatrasPayload
  PayloadFatrasExample.cpp)
target_link_libraries(
  ActsSimFatrasPayload
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorContextual)

# TGeo detector
add_executable(
  ActsSimFatrasTGeo
  TGeoFatrasExample.cpp)
target_link_libraries(
  ActsSimFatrasTGeo
  PRIVATE ActsExamplesFatrasCommon ActsExamplesDetectorTGeo)

install(
  TARGETS
    ActsSimFatrasGeneric
    ActsSimFatrasAligned
    ActsSimFatrasPayload
    ActsSimFatrasTGeo
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
