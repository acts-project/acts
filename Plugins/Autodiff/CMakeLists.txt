add_library(
  ActsPluginAutodiff INTERFACE)

target_include_directories(
  ActsPluginAutodiff
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)
 
# the target name is different when using find_package(...) or add_subdirectory(...)
if(ACTS_USE_SYSTEM_AUTODIFF)
    set(AUTODIFF_TARGET autodiff::autodiff)
else()
    set(AUTODIFF_TARGET autodiff)
endif()
    
target_link_libraries(
  ActsPluginAutodiff INTERFACE ActsCore ${AUTODIFF_TARGET})

install(
  TARGETS ActsPluginAutodiff
  EXPORT ActsPluginAutodiffTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
