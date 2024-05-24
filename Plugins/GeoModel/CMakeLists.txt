include(FetchContent)

add_library(
  ActsPluginGeoModel SHARED
  src/GeoModelConversionError.cpp
  src/GeoModelReader.cpp
  src/GeoModelDetectorElement.cpp
  src/GeoModelDetectorSurfaceFactory.cpp
  src/GeoModelSurfaceConverter.cpp
  src/converters/GeoBoxConverter.cpp
  src/converters/GeoTrdConverter.cpp
  src/converters/GeoTubeConverter.cpp
  src/converters/GeoIntersectionAnnulusConverter.cpp)
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
