add_library(
  ActsExamplesMaterialMapping SHARED
  src/MaterialMapping.cpp
  src/CoreMaterialMapping.cpp)
target_include_directories(
  ActsExamplesMaterialMapping
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesMaterialMapping
  PUBLIC ActsCore ActsExamplesFramework)

install(
  TARGETS ActsExamplesMaterialMapping
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
