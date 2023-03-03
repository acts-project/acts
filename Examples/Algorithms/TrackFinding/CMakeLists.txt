set(SOURCES 
  src/AmbiguityResolutionAlgorithm.cpp
  src/SeedingAlgorithm.cpp
  src/SeedingOrthogonalAlgorithm.cpp
  src/SpacePointMaker.cpp
  src/TrackFindingAlgorithm.cpp
  src/TrackFindingAlgorithmFunction.cpp
  src/HoughTransformSeeder.cpp
  src/TrackParamsEstimationAlgorithm.cpp
)
if(ACTS_BUILD_PLUGIN_ONNX)
  list(APPEND SOURCES 
    src/AmbiguityResolutionMLAlgorithm.cpp   
  )
endif()
add_library(
  ActsExamplesTrackFinding SHARED
  ${SOURCES}
)

target_include_directories(
  ActsExamplesTrackFinding
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTrackFinding
  PUBLIC
    ActsCore
    ActsExamplesFramework
    ActsExamplesIoJson
    ActsExamplesMagneticField)

if(ACTS_BUILD_PLUGIN_ONNX)
  target_compile_definitions(
    ActsExamplesTrackFinding PUBLIC ACTS_PLUGIN_ONNX)
  target_link_libraries(
    ActsExamplesTrackFinding PUBLIC ActsPluginOnnx)
endif()

install(
  TARGETS ActsExamplesTrackFinding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
