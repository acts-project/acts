set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesGenerators
  ActsExamplesGeneratorsPythia8
  ActsExamplesTruthTracking
  ActsExamplesVertexing
  ActsExamplesIoCsv
  ActsExamplesIoRoot
  Boost::program_options)

# Finder example
add_executable(ActsRecVertexFinder VertexFinderExample.cpp)
target_link_libraries(ActsRecVertexFinder PRIVATE ${_common_libraries})

# Fitter example
add_executable(ActsRecVertexFitter VertexFitterExample.cpp)
target_link_libraries(ActsRecVertexFitter PRIVATE ${_common_libraries})

# Writer example
add_executable(ActsRecVertexWriter VertexAndTracksWriterExample.cpp)
target_link_libraries(ActsRecVertexWriter PRIVATE ${_common_libraries})

# Reader example
add_executable(ActsRecVertexReader VertexAndTracksReaderAndFinderExample.cpp)
target_link_libraries(ActsRecVertexReader PRIVATE ${_common_libraries})

install(
  TARGETS
    ActsRecVertexFitter
    ActsRecVertexFinder
    ActsRecVertexWriter
    ActsRecVertexReader
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
