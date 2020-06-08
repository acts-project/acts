add_library(
  ActsTGeoPlugin SHARED
  src/TGeoDetectorElement.cpp
  src/TGeoLayerBuilder.cpp
  src/TGeoParser.cpp
  src/TGeoPrimitivesHelper.cpp
  src/TGeoSurfaceConverter.cpp)
target_include_directories(
  ActsTGeoPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsTGeoPlugin
  PUBLIC ActsCore ActsIdentificationPlugin ROOT::Geom)

install(
  TARGETS ActsTGeoPlugin
  EXPORT ActsTGeoPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
