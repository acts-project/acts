add_library(
  ActsExamplesDigitization SHARED
  src/DigitizationAlgorithm.cpp
  src/PlanarSteppingAlgorithm.cpp
  src/SmearingAlgorithm.cpp
  src/DigitizationOptions.cpp
  src/JsonDigitizationConfig.cpp)
target_include_directories(
  ActsExamplesDigitization
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDigitization
  PRIVATE
    ActsCore ActsPluginDigitization ActsPluginIdentification
    ActsExamplesFramework ActsExamplesIoJson
    Boost::program_options)

install(
  TARGETS ActsExamplesDigitization
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
