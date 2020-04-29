set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesIoObj
  ActsExamplesIoCsv
  ActsExamplesIoJson
  ActsExamplesIoRoot
  Boost::program_options)

# Empty detector
add_executable(
  ActsExampleGeometryEmpty
  EmptyGeometryExample.cpp)
target_link_libraries(
  ActsExampleGeometryEmpty
  PRIVATE ${_common_libraries} ActsExamplesDetectorEmpty)

# Generic detector
add_executable(
  ActsExampleGeometryGeneric
  GenericGeometryExample.cpp)
target_link_libraries(
  ActsExampleGeometryGeneric
  PRIVATE ${_common_libraries} ActsExamplesDetectorGeneric)

# Generic detector with IOV based alignment
add_executable(
  ActsExampleGeometryAligned
  AlignedGeometryExample.cpp)
target_link_libraries(
  ActsExampleGeometryAligned
  PRIVATE ${_common_libraries} ActsExamplesDetectorContextual)

# Generic detector with Payload based alignment
add_executable(
  ActsExampleGeometryPayload
  PayloadGeometryExample.cpp)
target_link_libraries(
  ActsExampleGeometryPayload
  PRIVATE ${_common_libraries} ActsExamplesDetectorContextual)

# TGeo detector
add_executable(
  ActsExampleGeometryTGeo
  TGeoGeometryExample.cpp)
target_link_libraries(
  ActsExampleGeometryTGeo
  PRIVATE ${_common_libraries} ActsExamplesDetectorTGeo)

install(
  TARGETS
    ActsExampleGeometryEmpty
    ActsExampleGeometryGeneric
    ActsExampleGeometryAligned
    ActsExampleGeometryPayload
    ActsExampleGeometryTGeo
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
