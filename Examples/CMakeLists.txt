add_component   (Framework ExamplesFramework)
add_subdirectory(Detectors)
add_subdirectory(Algorithms)
add_subdirectory(Io)
add_subdirectory_if(Scripts ACTS_BUILD_ANALYSIS_APPS)

add_subdirectory(HelloWorld)

propagate_components_to_parent()
