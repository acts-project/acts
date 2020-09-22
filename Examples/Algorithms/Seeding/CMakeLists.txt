add_library(
  ActsExamplesSeeding SHARED
  src/SeedingAlgorithm.cpp)

target_include_directories(
  ActsExamplesSeeding
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesSeeding
  PRIVATE
    ActsCore ActsPluginDigitization ActsPluginIdentification
    ActsExamplesFramework  ActsExamplesIoCsv
    Boost::program_options)

install(
  TARGETS ActsExamplesSeeding
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
