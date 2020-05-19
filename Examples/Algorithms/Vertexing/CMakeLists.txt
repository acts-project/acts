add_library(
  ActsExamplesVertexing SHARED
  src/VertexFitAlgorithm.cpp
  src/VertexFindingAlgorithm.cpp)
target_include_directories(
  ActsExamplesVertexing
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesVertexing
  PUBLIC ActsCore ActsExamplesFramework
  PRIVATE ActsExamplesTruthTracking)

install(
  TARGETS ActsExamplesVertexing
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
