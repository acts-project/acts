add_library(
  ActsPluginMlpack SHARED
  # header files
  include/Acts/Plugins/Mlpack/AmbiguityDBScanClustering.hpp)

target_include_directories(
  ActsPluginMlpack
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(
  ActsPluginMlpack
  PUBLIC ActsCore)

target_include_directories(ActsPluginMlpack 
  SYSTEM INTERFACE ${mlpack_INCLUDE_DIR})

set_target_properties(ActsPluginMlpack PROPERTIES LINKER_LANGUAGE CXX)

install(
  TARGETS ActsPluginMlpack
  EXPORT ActsPluginMlpackTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
