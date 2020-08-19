include(ExternalProject)

add_library(
  ActsAutodiffPlugin INTERFACE)

add_dependencies(
  ActsAutodiffPlugin autodiff)

target_include_directories(
  ActsAutodiffPlugin
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/thirdparty/autodiff>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    
target_link_libraries(
  ActsAutodiffPlugin INTERFACE ActsCore)
  
target_compile_definitions(
  ActsAutodiffPlugin INTERFACE AUTODIFF_PLUGIN_ENABLED)

install(
  TARGETS ActsAutodiffPlugin
  EXPORT ActsAutodiffPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
