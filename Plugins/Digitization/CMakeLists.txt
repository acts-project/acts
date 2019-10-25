add_library(
  ActsDigitizationPlugin SHARED
  src/CartesianSegmentation.cpp
  src/DigitizationModule.cpp
  src/PlanarModuleStepper.cpp)
target_include_directories(
  ActsDigitizationPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsDigitizationPlugin
  PUBLIC ActsCore)

install(
  TARGETS ActsDigitizationPlugin
  EXPORT ActsDigitizationPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
