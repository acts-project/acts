add_library(ActsIdentificationPlugin INTERFACE)
target_include_directories(
  ActsIdentificationPlugin
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(ActsIdentificationPlugin INTERFACE ActsCore)

install(
  TARGETS ActsIdentificationPlugin
  EXPORT ActsIdentificationPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
