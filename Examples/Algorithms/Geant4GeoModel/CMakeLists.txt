add_library(
  ActsExamplesGeant4GeoModel SHARED
  src/GdmlDetectorConstruction.cpp)
target_compile_definitions(
  ActsExamplesGeant4GeoModel
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsExamplesGeant4GeoModel
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(
  ActsExamplesGeant4GeoModel
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesGeant4GeoModel
  PUBLIC GeoModelCore::GeoModelKernel GeoModel2G4 ${Geant4_LIBRARIES}
  PRIVATE Boost::program_options)

install(
  TARGETS ActsExamplesGeant4GeoModel
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
