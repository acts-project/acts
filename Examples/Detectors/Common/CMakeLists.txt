add_library(
  ActsExamplesDetectorCommons SHARED
    src/Detector.cpp)
target_include_directories(
  ActsExamplesDetectorCommons
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorCommons
  PUBLIC
    ActsCore
    ActsExamplesFramework)

install(
  TARGETS ActsExamplesDetectorCommons
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
