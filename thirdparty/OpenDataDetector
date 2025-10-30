include(FetchContent)

if(NOT TARGET DD4hep::DDCore OR NOT TARGET DD4hep::DDDetectors)
    message(
        FATAL_ERROR
        "ACTS_BUILD_ODD requires DD4hep::DDCore and DD4hep::DDDetectors to be available as CMake targets."
    )
endif()

set(ACTS_ODD_GIT_REPOSITORY
    "https://gitlab.cern.ch/acts/OpenDataDetector.git"
    CACHE STRING
    "Git repository used to fetch OpenDataDetector"
)
set(ACTS_ODD_GIT_TAG
    "v4.0.3"
    CACHE STRING
    "Git tag or commit used for OpenDataDetector"
)

FetchContent_Declare(
    ODD
    GIT_REPOSITORY ${ACTS_ODD_GIT_REPOSITORY}
    GIT_TAG ${ACTS_ODD_GIT_TAG}
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(ODD)

FetchContent_GetProperties(
    ODD
    SOURCE_DIR ODD_SOURCE_DIR
    BINARY_DIR ODD_BINARY_DIR
)

set(ODD_FACTORY_DIR "${ODD_BINARY_DIR}/factory" CACHE INTERNAL "" FORCE)
set(ODD_SOURCE_DIR "${ODD_SOURCE_DIR}" CACHE INTERNAL "" FORCE)

message(STATUS "ODD source dir: ${ODD_SOURCE_DIR}")
message(STATUS "ODD binary dir: ${ODD_BINARY_DIR}")
message(STATUS "ODD factory dir: ${ODD_FACTORY_DIR}")
