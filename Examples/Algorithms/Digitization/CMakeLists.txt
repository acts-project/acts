add_library(
  ActsExamplesDigitization SHARED
  src/DigitizationAlgorithm.cpp
  src/DigitizationConfig.cpp
  src/MeasurementCreation.cpp
  src/PlanarSteppingAlgorithm.cpp
  src/ModuleClusters.cpp)
target_include_directories(
  ActsExamplesDigitization
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDigitization
  PUBLIC
    ActsCore ActsPluginIdentification ActsExamplesFramework)

install(
  TARGETS ActsExamplesDigitization
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
