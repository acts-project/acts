add_library(
  ActsExamplesGeant4DD4hep SHARED
  src/DD4hepDetectorConstruction.cpp)
target_compile_definitions(
  ActsExamplesGeant4DD4hep
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsExamplesGeant4DD4hep
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(
  ActsExamplesGeant4DD4hep
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesGeant4DD4hep
  PUBLIC 
  ActsExamplesGeant4 
  ${Geant4_LIBRARIES})

if(${DD4hep_VERSION} VERSION_LESS 1.11)
  target_include_directories(
    ActsExamplesGeant4DD4hep
    PRIVATE ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(
    ActsExamplesGeant4DD4hep
    PRIVATE ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDG4_LIBRARY})
else()
  target_link_libraries(
    ActsExamplesGeant4DD4hep
    PUBLIC DD4hep::DDCore DD4hep::DDG4)
endif()

install(
  TARGETS ActsExamplesGeant4DD4hep
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
