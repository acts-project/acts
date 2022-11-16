include(ActsTargetLinkLibrariesSystem)

add_library(
  ActsPluginJson SHARED
  src/MaterialMapJsonConverter.cpp
  src/AlgebraJsonConverter.cpp
  src/MaterialJsonConverter.cpp
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
  ActsPluginJson PUBLIC ActsCore)
acts_target_link_libraries_system(
  ActsPluginJson PUBLIC nlohmann_json::nlohmann_json)

install(
  TARGETS ActsPluginJson
  EXPORT ActsPluginJsonTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
