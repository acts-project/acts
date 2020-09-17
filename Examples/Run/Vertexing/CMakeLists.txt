set(_common_libraries
  ActsCore
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesGenerators
  ActsExamplesGeneratorsPythia8
  ActsExamplesTruthTracking
  ActsExamplesVertexing
  Boost::program_options)

add_executable(ActsExampleAdaptiveMultiVertexFinder AdaptiveMultiVertexFinderExample.cpp)
target_link_libraries(ActsExampleAdaptiveMultiVertexFinder PRIVATE ${_common_libraries})
add_executable(ActsExampleIterativeVertexFinder IterativeVertexFinderExample.cpp)
target_link_libraries(ActsExampleIterativeVertexFinder PRIVATE ${_common_libraries})
add_executable(ActsExampleVertexFitter VertexFitterExample.cpp)
target_link_libraries(ActsExampleVertexFitter PRIVATE ${_common_libraries})
add_executable(ActsTutorialVertexFinder TutorialVertexFinderExample.cpp)
target_link_libraries(ActsTutorialVertexFinder PRIVATE ${_common_libraries})

# Writer example
# add_executable(ActsExampleVertexWriter VertexAndTracksWriterExample.cpp)
# target_link_libraries(ActsExampleVertexWriter PRIVATE ${_common_libraries})
#
# # Reader example
# add_executable(ActsExampleVertexReader VertexAndTracksReaderAndFinderExample.cpp)
# target_link_libraries(ActsExampleVertexReader PRIVATE ${_common_libraries})

# Tutorial AMVF example

install(
  TARGETS
    ActsExampleAdaptiveMultiVertexFinder
    ActsExampleIterativeVertexFinder
    ActsExampleVertexFitter
    ActsTutorialVertexFinder
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
