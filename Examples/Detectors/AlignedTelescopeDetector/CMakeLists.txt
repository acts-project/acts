add_library(
  ActsExamplesDetectorAlignedTelescope SHARED
  src/AlignedTelescopeDetector.cpp
  src/AlignedTelescopeDetectorElement.cpp
  src/BuildAlignedTelescopeDetector.cpp)
target_include_directories(
  ActsExamplesDetectorAlignedTelescope
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorAlignedTelescope
  PUBLIC ActsCore ActsExamplesFramework)

install(
  TARGETS ActsExamplesDetectorAlignedTelescope
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
