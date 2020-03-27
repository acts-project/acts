add_library(
  ActsExamplesDetectorEmpty SHARED
  src/EmptyDetector.cpp)
target_include_directories(
  ActsExamplesDetectorEmpty
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorEmpty
  PUBLIC ActsCore ActsExamplesFramework ActsExamplesDetectorsCommon)

install(
  TARGETS ActsExamplesDetectorEmpty
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
