add_library(
  ActsExamplesIoJson SHARED
  src/JsonDigitizationConfig.cpp
  src/JsonMaterialWriter.cpp
  src/JsonSurfacesWriter.cpp)
target_include_directories(
  ActsExamplesIoJson
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoJson
  PUBLIC ActsCore ActsPluginJson ActsExamplesFramework ActsExamplesDigitization)

install(
  TARGETS ActsExamplesIoJson
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
