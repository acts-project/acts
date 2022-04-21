add_library(
  ActsPluginIdentification INTERFACE)
target_include_directories(
  ActsPluginIdentification
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(ActsPluginIdentification INTERFACE ActsCore)

install(
  TARGETS ActsPluginIdentification
  EXPORT ActsPluginIdentificationTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
