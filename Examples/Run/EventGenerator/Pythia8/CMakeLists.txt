add_executable(
  ActsExamplePythia8
  Pythia8.cpp)
target_link_libraries(
  ActsExamplePythia8
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesGenerators ActsExamplesGeneratorsPythia8
    ActsExamplesIoCsv ActsExamplesIoRoot ActsExamplesPrinters
    Boost::program_options)

install(
  TARGETS ActsExamplePythia8
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
