add_library(
  ActsPluginAutodiff INTERFACE)

target_include_directories(
  ActsPluginAutodiff
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)
    
target_link_libraries(
  ActsPluginAutodiff INTERFACE ActsCore autodiff::autodiff)

install(
  TARGETS ActsPluginAutodiff
  EXPORT ActsPluginAutodiffTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
