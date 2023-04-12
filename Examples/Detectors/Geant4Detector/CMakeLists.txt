add_library(
  ActsExamplesDetectorGeant4 SHARED
  src/Geant4Detector.cpp
  src/MockupSectorBuilder.cpp)
target_include_directories(
    ActsExamplesDetectorGeant4
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
    ActsExamplesDetectorGeant4
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesGeant4 ActsPluginGeant4)

add_executable(
  ActsExampleTestMockupBuilder
  src/TestMockupBuilder.cpp)
target_include_directories(
  ActsExampleTestMockupBuilder
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExampleTestMockupBuilder
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesGeant4 ActsPluginGeant4
  ActsExamplesFramework
  ActsExamplesGeant4
  ActsExamplesIoRoot
  ActsExamplesDetectorGeant4
)

install(
  TARGETS ActsExamplesDetectorGeant4 ActsExampleTestMockupBuilder
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
