add_library(
  ActsPluginJson SHARED
  src/AlgebraJsonConverter.cpp
  src/DetectorJsonConverter.cpp
  src/DetrayJsonHelper.cpp
  src/GridJsonConverter.cpp
  src/MaterialMapJsonConverter.cpp
  src/MaterialJsonConverter.cpp
  src/PortalJsonConverter.cpp
  src/SurfaceBoundsJsonConverter.cpp
  src/SurfaceJsonConverter.cpp
  src/UtilitiesJsonConverter.cpp
  src/VolumeBoundsJsonConverter.cpp
  src/VolumeJsonConverter.cpp)
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
