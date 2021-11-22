# Json based
add_executable(
  ActsExampleTracerDetector
  TracerDetectorExample.cpp)
target_link_libraries(
  ActsExampleTracerDetector
  PRIVATE ActsPluginJson ActsExamplesFramework ActsExamplesPropagation)


# Json based
add_executable(
  ActsExampleTracerPropagation
  TracerPropagationExample.cpp)
target_link_libraries(
  ActsExampleTracerPropagation
  PRIVATE ActsExamplesFramework ActsExamplesPropagation)