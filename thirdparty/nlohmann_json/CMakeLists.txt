cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building nlohmann_json as part of the ACTS project" )


# Declare where to get VecMem from.
set( ACTS_NLOHMANN_JSON_GIT_REPOSITORY "https://github.com/nlohmann/json.git"
   CACHE STRING "Git repository to take nlohmann_json from" )
set( ACTS_NLOHMANN_JSON_GIT_TAG "v3.10.5" CACHE STRING "Version of nlohmann_json to build" )
mark_as_advanced( ACTS_NLOHMANN_JSON_GIT_REPOSITORY ACTS_NLOHMANN_JSON_GIT_TAG )
FetchContent_Declare( nlohmann_json
   GIT_REPOSITORY "${ACTS_NLOHMANN_JSON_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_NLOHMANN_JSON_GIT_TAG}")

# Now set up its build.
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install ON CACHE INTERNAL "")
FetchContent_MakeAvailable( nlohmann_json )
