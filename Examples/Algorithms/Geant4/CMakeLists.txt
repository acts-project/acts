add_library(
  ActsExamplesGeant4 SHARED
  src/EventAction.cpp
  src/GeantinoRecording.cpp
  src/Geant4Options.cpp
  src/PrimaryGeneratorAction.cpp
  src/RunAction.cpp
  src/SteppingAction.cpp)
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
  PUBLIC ActsCore ActsExamplesFramework ${Geant4_LIBRARIES}
  PRIVATE Boost::program_options)

install(
  TARGETS ActsExamplesGeant4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

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
