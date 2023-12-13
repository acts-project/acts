cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building nlohmann_json as part of the ACTS project" )


# Declare where to get VecMem from.
set( ACTS_NLOHMANN_JSON_GIT_TAG "v3.10.5" CACHE STRING "Version of nlohmann_json to build" )
set( ACTS_NLOHMANN_JSON_SHA_256 "d6c65aca6b1ed68e7a182f4757257b107ae403032760ed6ef121c9d55e81757d" CACHE STRING "SHA256 hash of the downloaded zip" )
mark_as_advanced( ACTS_NLOHMANN_JSON_GIT_REPOSITORY ACTS_NLOHMANN_JSON_GIT_TAG )
FetchContent_Declare( nlohmann_json
   URL "https://github.com/nlohmann/json/releases/download/${ACTS_NLOHMANN_JSON_GIT_TAG}/json.tar.xz"
   URL_HASH SHA256=${ACTS_NLOHMANN_JSON_SHA_256})

# Now set up its build.
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install ON CACHE INTERNAL "")
FetchContent_MakeAvailable( nlohmann_json )
