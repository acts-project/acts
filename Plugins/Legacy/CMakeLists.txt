add_library(
  ActsLegacy INTERFACE)
target_compile_features(
  ActsLegacy
  INTERFACE cxx_std_17)
target_include_directories(
  ActsLegacy
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsLegacy
  INTERFACE Boost::boost Eigen3::Eigen)

install(
  TARGETS ActsLegacy
  EXPORT ActsLegacyTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
