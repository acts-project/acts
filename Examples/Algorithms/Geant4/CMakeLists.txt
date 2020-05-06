add_library(
  ActsExamplesGeant4 SHARED
  src/DD4hepToG4Svc.cpp
  src/GeantinoRecording.cpp
  src/GeoConstruction.cpp
  src/MMDetectorConstruction.cpp
  src/MMEventAction.cpp
  src/MMPrimaryGeneratorAction.cpp
  src/MMRunAction.cpp
  src/MMSteppingAction.cpp)
target_compile_definitions(
  ActsExamplesGeant4
  PUBLIC ${Geant4_DEFINITIONS})
target_include_directories(
  ActsExamplesGeant4
  SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
target_include_directories(
  ActsExamplesGeant4
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesGeant4
  PUBLIC
    ActsCore ActsExamplesDetectorDD4hep ActsExamplesFramework
    ${Geant4_LIBRARIES})

if(${DD4hep_VERSION} VERSION_LESS 1.11)
  target_include_directories(
    ActsExamplesGeant4
    PUBLIC ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(
    ActsExamplesGeant4
    PUBLIC ${DD4hep_DDCORE_LIBRARY} ${DD4hep_DDG4_LIBRARY})
else()
  target_link_libraries(
    ActsExamplesGeant4
    PUBLIC DD4hep::DDCore DD4hep::DDG4)
endif()

install(
  TARGETS ActsExamplesGeant4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
