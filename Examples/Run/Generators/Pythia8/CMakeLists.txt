add_library(ActsExamplesGeneratorsPythia8Options SHARED
  Options/src/Pythia8Options.cpp)
target_include_directories(
  ActsExamplesGeneratorsPythia8Options
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Options/include>)
target_link_libraries(
  ActsExamplesGeneratorsPythia8Options
  PRIVATE
    ActsExamplesFramework
    ActsExamplesGenerators
    ActsExamplesGeneratorsPythia8
    Boost::program_options)

add_executable(
  ActsExamplePythia8
  Pythia8.cpp)
target_link_libraries(
  ActsExamplePythia8
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesGenerators ActsExamplesGeneratorsPythia8
    ActsExamplesIoCsv ActsExamplesIoRoot ActsExamplesPrinters
    ActsExamplesGeneratorsPythia8Options)

install(
  TARGETS ActsExamplePythia8
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
