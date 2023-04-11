add_library(
  ActsPluginmlpack SHARED
  # header files
  include/Acts/Plugins/mlpack/AmbiguityDBScanClustering.hpp)

target_include_directories(
  ActsPluginmlpack
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(
  ActsPluginmlpack
  PUBLIC ActsCore)

set_target_properties(ActsPluginmlpack PROPERTIES LINKER_LANGUAGE CXX)

install(
  TARGETS ActsPluginmlpack
  EXPORT ActsPluginmlpackTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
