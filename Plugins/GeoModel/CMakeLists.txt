include(FetchContent)
add_library(
  ActsPluginGeoModel SHARED
  src/GeoModelToDetVol.cpp)
target_include_directories(
    ActsPluginGeoModel
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(
      ActsPluginGeoModel
      SYSTEM PUBLIC)    
target_link_libraries(
    ActsPluginGeoModel
  PUBLIC ActsCore GeoModelCore::GeoModelKernel)
install(
  TARGETS ActsPluginGeoModel
  EXPORT ActsPluginGeoModelTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})