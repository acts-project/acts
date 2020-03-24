add_library(
  ActsDD4hepPlugin SHARED
  src/ActsExtension.cpp
  src/ConvertDD4hepDetector.cpp
  src/ConvertDD4hepMaterial.cpp
  src/DD4hepDetectorElement.cpp
  src/DD4hepLayerBuilder.cpp
  src/DD4hepVolumeBuilder.cpp)
target_include_directories(
  ActsDD4hepPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    ${DD4hep_INCLUDE_DIRS})
target_link_libraries(
  ActsDD4hepPlugin
  PUBLIC ActsCore ActsTGeoPlugin)

if(${DD4hep_VERSION} VERSION_LESS 1.11)
  target_link_libraries(ActsDD4hepPlugin PUBLIC ${DD4hep_DDCORE_LIBRARY})
else()
  target_link_libraries(ActsDD4hepPlugin PUBLIC DD4hep::DDCore)
endif()

install(
  TARGETS ActsDD4hepPlugin
  EXPORT ActsDD4hepPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
