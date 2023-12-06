add_library(
  ActsExamplesDetectorF2Telescope SHARED
  src/F2TelescopeDetector.cpp
  src/F2TelescopeDetectorElement.cpp
  src/BuildF2TelescopeDetector.cpp)
target_include_directories(
  ActsExamplesDetectorF2Telescope
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorF2Telescope
  PUBLIC ActsCore ActsExamplesFramework)

install(
  TARGETS ActsExamplesDetectorF2Telescope
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
