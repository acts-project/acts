cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building ActsDD4hep glue library as part of the ACTS project" )


# Declare where to get VecMem from.
set( ACTS_ACTSDD4HEP_GIT_REPOSITORY "https://github.com/acts-project/acts-dd4hep.git"
   CACHE STRING "Git repository to take ActsDD4hep glue library from from" )
set( ACTS_ACTSDD4HEP_GIT_TAG "v1.0.0" CACHE STRING "Version of ActsDD4hep glue library to build" )
mark_as_advanced( ACTS_ACTSDD4HEP_GIT_REPOSITORY ACTS_ACTSDD4HEP_GIT_TAG )
FetchContent_Declare( actsdd4hep
   GIT_REPOSITORY "${ACTS_ACTSDD4HEP_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_ACTSDD4HEP_GIT_TAG}" )

FetchContent_MakeAvailable( actsdd4hep )
