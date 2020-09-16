add_library(
  ActsPluginDigitization SHARED
  src/CartesianSegmentation.cpp
  src/DigitizationModule.cpp
  src/PlanarModuleStepper.cpp)
target_include_directories(
  ActsPluginDigitization
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginDigitization
  PUBLIC ActsCore)

install(
  TARGETS ActsPluginDigitization
  EXPORT ActsPluginDigitizationTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
