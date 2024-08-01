cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building dfelibs as part of the ACTS project" )

set( DFELIBS_VERSION "v${_acts_dfelibs_version}")

# Declare where to get dfelibs from.
FetchContent_Declare( dfelibs ${ACTS_DFELIBS_SOURCE} )

# Options used in the build of dfelibs.
set( dfelibs_BUILD_EXAMPLES OFF CACHE BOOL
   "Turn off the build of the dfelibs examples" )
set( dfelibs_BUILD_UNITTESTS OFF CACHE BOOL
   "Turn off the build of the dfelibs unit tests" )
set( dfelibs_ENABLE_INSTALL ON CACHE BOOL
   "Have dfelibs be installed together with the rest of the project" )

# Now set up its build.
FetchContent_MakeAvailable( dfelibs )
