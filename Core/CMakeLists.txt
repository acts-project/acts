configure_file(
    ActsVersion.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/Acts/ActsVersion.hpp
)

# source files will be added later
acts_add_library(Core "")

target_sources(ActsCore PRIVATE src/ActsVersion.cpp)

target_compile_features(ActsCore PUBLIC ${ACTS_CXX_STANDARD_FEATURE})

target_include_directories(
    ActsCore
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        # for the generated version header
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(ActsCore PUBLIC Boost::boost Eigen3::Eigen)

if(ACTS_PARAMETER_DEFINITIONS_HEADER)
    target_compile_definitions(
        ActsCore
        PUBLIC
            -DACTS_PARAMETER_DEFINITIONS_HEADER="${ACTS_PARAMETER_DEFINITIONS_HEADER}"
    )
endif()

if(ACTS_SOURCELINK_SBO_SIZE)
    target_compile_definitions(
        ActsCore
        PUBLIC -DACTS_SOURCELINK_SBO_SIZE=${ACTS_SOURCELINK_SBO_SIZE}
    )
endif()

if(ACTS_LOG_FAILURE_THRESHOLD)
    message(
        STATUS
        "Enable log failure threshold, set to ${ACTS_LOG_FAILURE_THRESHOLD}"
    )
    target_compile_definitions(
        ActsCore
        PUBLIC
            -DACTS_LOG_FAILURE_THRESHOLD=${ACTS_LOG_FAILURE_THRESHOLD}
            -DACTS_ENABLE_LOG_FAILURE_THRESHOLD
    )
else()
    if(ACTS_ENABLE_LOG_FAILURE_THRESHOLD)
        message(STATUS "Enable log failure threshold")
        target_compile_definitions(
            ActsCore
            PUBLIC -DACTS_ENABLE_LOG_FAILURE_THRESHOLD
        )
    endif()
endif()

if(ACTS_ENABLE_CPU_PROFILING)
    message(STATUS "added lprofiler")

    if(NOT ACTS_GPERF_INSTALL_DIR STREQUAL "")
        find_library(PROFILER_LIBRARY profiler HINTS ${ACTS_GPERF_INSTALL_DIR})

        target_link_libraries(ActsCore PUBLIC ${PROFILER_LIBRARY})
    else()
        target_link_libraries(ActsCore PUBLIC -lprofiler)
    endif()

    target_link_options(ActsCore PUBLIC "LINKER:-no-as-needed")
endif()

if(ACTS_ENABLE_MEMORY_PROFILING)
    message(STATUS "added ltcmalloc")

    if(NOT ACTS_GPERF_INSTALL_DIR STREQUAL "")
        find_library(TCMALLOC_LIBRARY tcmalloc HINTS ${ACTS_GPERF_INSTALL_DIR})

        target_link_libraries(ActsCore PUBLIC ${TCMALLOC_LIBRARY})
    else()
        target_link_libraries(ActsCore PUBLIC -ltcmalloc)
    endif()
endif()

install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install version header based on git repository status
if(_acts_is_git_repo)
    # For git repositories, generate header with hash at install time
    install(
        CODE
            "
        # Get git hash at install time
        execute_process(
            COMMAND git rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE _install_git_hash
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND git rev-parse --short HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE _install_git_hash_short
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        execute_process(
            COMMAND git diff-index --quiet HEAD --
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            RESULT_VARIABLE _git_is_dirty
            ERROR_QUIET
        )

        if(NOT _git_is_dirty EQUAL 0)
            set(_install_git_hash \"\${_install_git_hash}-dirty\")
            set(_install_git_hash_short \"\${_install_git_hash_short}-dirty\")
        endif()

        # Configure the version file with hash and version numbers
        set(_acts_commit_hash \"\\\"\${_install_git_hash}\\\"\")
        set(_acts_commit_hash_short \"\\\"\${_install_git_hash_short}\\\"\")
        set(PROJECT_VERSION_MAJOR \"${PROJECT_VERSION_MAJOR}\")
        set(PROJECT_VERSION_MINOR \"${PROJECT_VERSION_MINOR}\")
        set(PROJECT_VERSION_PATCH \"${PROJECT_VERSION_PATCH}\")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/ActsVersion.hpp.in
        ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/Acts/ActsVersion.hpp
        @ONLY
    )
"
    )
else()
    # For non-git sources (e.g., tarball), just copy the build-time header
    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/include/Acts/ActsVersion.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Acts
    )
endif()

# target source files are added separately
add_subdirectory(src/EventData)
add_subdirectory(src/Definitions)
add_subdirectory(src/Geometry)
add_subdirectory(src/MagneticField)
add_subdirectory(src/Material)
add_subdirectory(src/Navigation)
add_subdirectory(src/Propagator)
add_subdirectory(src/Seeding)
add_subdirectory(src/Seeding2)
add_subdirectory(src/SpacePointFormation2)
add_subdirectory(src/Surfaces)
add_subdirectory(src/TrackFinding)
add_subdirectory(src/TrackFitting)
add_subdirectory(src/Utilities)
add_subdirectory(src/Vertexing)
add_subdirectory(src/Visualization)
add_subdirectory(src/AmbiguityResolution)

acts_compile_headers(Core GLOB include/Acts/**/*.hpp)
