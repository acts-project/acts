add_library(
  ActsExamplesGeant4 SHARED
  src/GdmlDetectorConstruction.cpp
  src/GeantinoRecording.cpp
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

install(
  TARGETS ActsExamplesGeant4
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
