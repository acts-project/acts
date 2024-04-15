cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building dfelibs as part of the ACTS project" )

set( DFELIBS_VERSION "v${_acts_dfelibs_version}")

# Declare where to get VecMem from.
set( ACTS_DFELIBS_GIT_REPOSITORY "https://github.com/acts-project/dfelibs.git"
   CACHE STRING "Git repository to take dfelibs from" )
set( ACTS_DFELIBS_GIT_TAG "${DFELIBS_VERSION}" CACHE STRING "Version of dfelibs to build" )
mark_as_advanced( ACTS_DFELIBS_GIT_REPOSITORY ACTS_DFELIBS_GIT_TAG )
FetchContent_Declare( dfelibs
   GIT_REPOSITORY "${ACTS_DFELIBS_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_DFELIBS_GIT_TAG}" )

# Options used in the build of dfelibs.
set( dfelibs_BUILD_EXAMPLES OFF CACHE BOOL
   "Turn off the build of the dfelibs examples" )
set( dfelibs_BUILD_UNITTESTS OFF CACHE BOOL
   "Turn off the build of the dfelibs unit tests" )
set( dfelibs_ENABLE_INSTALL ON CACHE BOOL
   "Have dfelibs be installed together with the rest of the project" )

# Now set up its build.
FetchContent_MakeAvailable( dfelibs )
