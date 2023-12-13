cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building nlohmann_json as part of the ACTS project" )


# Declare where to get VecMem from.
set( ACTS_NLOHMANN_JSON_GIT_TAG "v3.10.5" CACHE STRING "Version of nlohmann_json to build" )
set( ACTS_NLOHMANN_JSON_SHA_256 "344be97b757a36c5b180f1c8162f6c5f6ebd760b117f6e64b77866e97b217280" CACHE STRING "SHA256 hash of the downloaded zip" )
mark_as_advanced( ACTS_NLOHMANN_JSON_GIT_REPOSITORY ACTS_NLOHMANN_JSON_GIT_TAG )
FetchContent_Declare( nlohmann_json
   URL "https://github.com/nlohmann/json/releases/download/${ACTS_NLOHMANN_JSON_GIT_TAG}/json.tar.xz"
   URL_HASH SHA256=${ACTS_NLOHMANN_JSON_SHA_256})

# Now set up its build.
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install ON CACHE INTERNAL "")
FetchContent_MakeAvailable( nlohmann_json )
