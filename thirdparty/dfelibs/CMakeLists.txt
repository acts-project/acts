cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building dfelibs as part of the ACTS project" )


# Declare where to get VecMem from.
set( ACTS_DFELIBS_GIT_REPOSITORY "https://github.com/acts-project/dfelibs.git"
   CACHE STRING "Git repository to take dfelibs from" )
set( ACTS_DFELIBS_GIT_TAG "v20200416" CACHE STRING "Version of dfelibs to build" )
mark_as_advanced( ACTS_DFELIBS_GIT_REPOSITORY ACTS_DFELIBS_GIT_TAG )
FetchContent_Declare( dfelibs
   GIT_REPOSITORY "${ACTS_DFELIBS_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_DFELIBS_GIT_TAG}" )

# Now set up its build.
FetchContent_MakeAvailable( dfelibs )
