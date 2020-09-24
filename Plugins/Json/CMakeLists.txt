add_library(
  ActsJsonPlugin SHARED
  src/JsonGeometryConverter.cpp)
target_include_directories(
  ActsJsonPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  ActsJsonPlugin PUBLIC ActsCore nlohmann_json::nlohmann_json)

install(
  TARGETS ActsJsonPlugin
  EXPORT ActsJsonPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
