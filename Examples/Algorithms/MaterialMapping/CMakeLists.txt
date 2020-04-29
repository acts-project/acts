add_library(
  ActsExamplesMaterialMapping SHARED
  src/MaterialMapping.cpp)
target_include_directories(
  ActsExamplesMaterialMapping
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesMaterialMapping
  PUBLIC ActsCore ActsExamplesFramework)

if(ACTS_BUILD_EXAMPLES_DD4HEP AND ACTS_BUILD_EXAMPLES_GEANT4)
  include(${Geant4_USE_FILE})
  target_sources(ActsExamplesMaterialMapping PRIVATE src/GeantinoRecording.cpp)
  target_link_libraries(
    ActsExamplesMaterialMapping
    PUBLIC ActsExamplesGeant4 ActsExamplesDD4hepG4)
endif()

install(
  TARGETS ActsExamplesMaterialMapping
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
