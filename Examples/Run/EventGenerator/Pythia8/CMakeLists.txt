add_executable(
  ActsGenPythia8
  Pythia8.cpp)
target_link_libraries(
  ActsGenPythia8
  PRIVATE
    ActsExamplesFramework ActsExamplesCommon
    ActsExamplesGenerators ActsExamplesGeneratorsPythia8
    ActsExamplesIoRoot ActsExamplesIoCsv ActsExamplesPrinters
    Boost::program_options)

install(
  TARGETS ActsGenPythia8
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
