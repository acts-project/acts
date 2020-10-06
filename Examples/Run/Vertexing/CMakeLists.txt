set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesTruthTracking
  ActsExamplesVertexing
  Boost::program_options)

add_executable(
  ActsExampleVertexFinderReader
  VertexFinderReaderExample.cpp)
target_link_libraries(
  ActsExampleVertexFinderReader
  PRIVATE ${_common_libraries} ActsExamplesIoCsv ActsExamplesPrinters)

install(
  TARGETS
    ActsExampleVertexFinderReader
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(Pythia8 ACTS_BUILD_EXAMPLES_PYTHIA8)
