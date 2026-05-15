add_subdirectory(MaterialMapping)
add_subdirectory(TrackingPerformance)

if(ACTS_BUILD_EXAMPLES_ROOT)
    add_executable(compareRootFiles compareRootFiles.C)
    target_link_libraries(
        compareRootFiles
        PUBLIC ROOT::Core ROOT::Tree ROOT::TreePlayer
    )
endif()
