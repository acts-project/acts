include(FetchContent)

add_library(
  ActsPluginGeoModel SHARED
  src/GeoModelReader.cpp)
target_include_directories(
    ActsPluginGeoModel
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
    ActsPluginGeoModel
  PUBLIC ActsCore GeoModelKernel GeoModelIO::GeoModelDBManager GeoModelIO::GeoModelRead)

install(
  TARGETS ActsPluginGeoModel
  EXPORT ActsPluginGeoModelTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
