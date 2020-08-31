add_library(
  ActsPluginJson SHARED
  src/JsonGeometryConverter.cpp)
target_include_directories(
  ActsPluginJson
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  ActsPluginJson PUBLIC ActsCore nlohmann_json::nlohmann_json)

install(
  TARGETS ActsPluginJson
  EXPORT ActsPluginJsonTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
