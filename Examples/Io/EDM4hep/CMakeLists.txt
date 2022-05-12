add_library(
  ActsExamplesIoEDM4hep SHARED
  src/test.cpp)
target_include_directories(
  ActsExamplesIoEDM4hep
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoEDM4hep
  PUBLIC ActsCore ActsPluginEDM4hep)

install(
  TARGETS ActsExamplesIoEDM4hep
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
