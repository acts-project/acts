add_library(
  ActsExamplesMultiWire SHARED
  src/MultiWireStructureBuilder.cpp
  src/MultiWireHelper.cpp)
target_include_directories(
    ActsExamplesMultiWire
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
    ActsExamplesMultiWire
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesGeant4 ActsExamplesDetectorGeant4 ActsPluginGeant4)

install(
  TARGETS ActsExamplesMultiWire
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
