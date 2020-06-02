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
add_executable(ActsRecIterativeVertexFinder IterativeVertexFinderExample.cpp)
target_link_libraries(ActsRecIterativeVertexFinder PRIVATE ${_common_libraries})

# Fitter example
add_executable(ActsRecVertexFitter VertexFitterExample.cpp)
target_link_libraries(ActsRecVertexFitter PRIVATE ${_common_libraries})

# Writer example
add_executable(ActsRecVertexWriter VertexAndTracksWriterExample.cpp)
target_link_libraries(ActsRecVertexWriter PRIVATE ${_common_libraries})

# Reader example
add_executable(ActsRecVertexReader VertexAndTracksReaderAndFinderExample.cpp)
target_link_libraries(ActsRecVertexReader PRIVATE ${_common_libraries})

# Tutorial AMVF example
add_executable(ActsTutorialAMVF TutorialAMVFExample.cpp)
target_link_libraries(ActsTutorialAMVF PRIVATE ${_common_libraries})

# Full AMVF example
add_executable(ActsRecAdaptiveMultiVertexFinder AdaptiveMultiVertexFinderExample.cpp)
target_link_libraries(ActsRecAdaptiveMultiVertexFinder PRIVATE ${_common_libraries})

install(
  TARGETS
    ActsRecVertexFitter
    ActsRecIterativeVertexFinder
    ActsRecVertexWriter
    ActsRecVertexReader
    ActsRecAdaptiveMultiVertexFinder
    ActsTutorialAMVF
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
