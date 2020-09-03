add_library(
  ActsExamplesIoJson SHARED
  src/JsonMaterialWriter.cpp)
target_include_directories(
  ActsExamplesIoJson
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoJson
  PUBLIC ActsCore ActsJsonPlugin ActsExamplesFramework)

install(
  TARGETS ActsExamplesIoJson
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
