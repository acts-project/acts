include(${Geant4_USE_FILE})

add_library(
  ActsExamplesGeant4 SHARED
  src/MMDetectorConstruction.cpp
  src/MMEventAction.cpp
  src/MMPrimaryGeneratorAction.cpp
  src/MMRunAction.cpp
  src/MMSteppingAction.cpp)
target_include_directories(
  ActsExamplesGeant4
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${Geant4_INCLUDE_DIRS})
target_link_libraries(
  ActsExamplesGeant4
  PUBLIC ActsCore ActsExamplesFramework ${Geant4_LIBRARIES}
  PRIVATE ROOT::Geom)

install(
  TARGETS ActsExamplesGeant4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
