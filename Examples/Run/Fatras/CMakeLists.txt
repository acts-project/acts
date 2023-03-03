# Common library that implements the core functionality
add_subdirectory(Common)

# Generic detector
add_executable(
  ActsExampleFatrasGeneric
  GenericFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasGeneric
  PRIVATE ActsExamplesFatrasCommon ActsExamplesCommon ActsExamplesIoJson)

# Generic detector with alignment
add_executable(
  ActsExampleFatrasAligned
  AlignedFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasAligned
  PRIVATE ActsExamplesFatrasCommon ActsExamplesCommon ActsExamplesIoJson)

# TGeo detector
add_executable(
  ActsExampleFatrasTGeo
  TGeoFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasTGeo
  PRIVATE ActsExamplesFatrasCommon ActsExamplesCommon ActsExamplesIoJson)

# Telescope detector
add_executable(
  ActsExampleFatrasTelescope
  TelescopeFatrasExample.cpp)
target_link_libraries(
  ActsExampleFatrasTelescope
  PRIVATE ActsExamplesFatrasCommon ActsExamplesCommon ActsExamplesIoJson)

install(
  TARGETS
    ActsExampleFatrasGeneric
    ActsExampleFatrasAligned
    ActsExampleFatrasTGeo
    ActsExampleFatrasTelescope
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
