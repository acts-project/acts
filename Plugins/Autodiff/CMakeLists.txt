add_library(
  ActsPluginAutodiff INTERFACE)

add_dependencies(
  ActsPluginAutodiff autodiff)

target_include_directories(
  ActsPluginAutodiff
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/thirdparty/autodiff>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    
target_link_libraries(
  ActsPluginAutodiff INTERFACE ActsCore)

install(
  TARGETS ActsPluginAutodiff
  EXPORT ActsPluginAutodiffTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
