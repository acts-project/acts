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
add_executable(ActsExampleIterativeVertexFinder IterativeVertexFinderExample.cpp)
target_link_libraries(ActsExampleIterativeVertexFinder PRIVATE ${_common_libraries})

# Fitter example
add_executable(ActsExampleVertexFitter VertexFitterExample.cpp)
target_link_libraries(ActsExampleVertexFitter PRIVATE ${_common_libraries})

# Writer example
add_executable(ActsExampleVertexWriter VertexAndTracksWriterExample.cpp)
target_link_libraries(ActsExampleVertexWriter PRIVATE ${_common_libraries})

# Reader example
add_executable(ActsExampleVertexReader VertexAndTracksReaderAndFinderExample.cpp)
target_link_libraries(ActsExampleVertexReader PRIVATE ${_common_libraries})

# Tutorial AMVF example
add_executable(ActsTutorialAMVF TutorialAMVFExample.cpp)
target_link_libraries(ActsTutorialAMVF PRIVATE ${_common_libraries})

# Full AMVF example
add_executable(ActsExampleAdaptiveMultiVertexFinder AdaptiveMultiVertexFinderExample.cpp)
target_link_libraries(ActsExampleAdaptiveMultiVertexFinder PRIVATE ${_common_libraries})

install(
  TARGETS
    ActsExampleVertexFitter
    ActsExampleIterativeVertexFinder
    ActsExampleVertexWriter
    ActsExampleVertexReader
    ActsExampleAdaptiveMultiVertexFinder
    ActsTutorialAMVF
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
