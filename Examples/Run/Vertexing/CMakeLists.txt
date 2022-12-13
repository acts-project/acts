set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesTruthTracking
  ActsExamplesVertexing)

add_executable(
  ActsExampleVertexFinderParticleReader
  ParticleReaderVertexingExample.cpp)
target_link_libraries(
  ActsExampleVertexFinderParticleReader
  PRIVATE ${_common_libraries} ActsExamplesIoCsv ActsExamplesPrinters)

  add_executable(
  ActsExampleVertexFinderTrackReaderPerformanceWriter
  TrackReaderVertexingPerformanceWriterExample.cpp)
target_link_libraries(
  ActsExampleVertexFinderTrackReaderPerformanceWriter
  PRIVATE ${_common_libraries} ActsExamplesIoRoot ActsExamplesPrinters)

install(
  TARGETS
    ActsExampleVertexFinderParticleReader
    ActsExampleVertexFinderTrackReaderPerformanceWriter
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(Pythia8 ACTS_BUILD_EXAMPLES_PYTHIA8)
