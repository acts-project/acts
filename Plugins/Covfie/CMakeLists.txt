add_library(
  ActsPluginCovfie INTERFACE)

target_include_directories(
  ActsPluginCovfie
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginCovfie
  INTERFACE
    ActsCore
    covfie::core)

install(
  TARGETS ActsPluginCovfie
  EXPORT ActsPluginCovfieTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
