cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Tell the user what's happening.
message( STATUS "Building autodiff as part of the ACTS project" )


# Declare where to get VecMem from.
set( ACTS_AUTODIFF_GIT_REPOSITORY "https://github.com/autodiff/autodiff.git"
   CACHE STRING "Git repository to take autodiff from" )
set( ACTS_AUTODIFF_GIT_TAG "v0.6.4" CACHE STRING "Version of autodiff to build" )
mark_as_advanced( ACTS_AUTODIFF_GIT_REPOSITORY ACTS_AUTODIFF_GIT_TAG )
FetchContent_Declare( autodiff
   GIT_REPOSITORY "${ACTS_AUTODIFF_GIT_REPOSITORY}"
   GIT_TAG "${ACTS_AUTODIFF_GIT_TAG}" )

# Now set up its build.
# This policy allows to set the options outside without changing autodiff's CMakeLists.txt
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(AUTODIFF_BUILD_TESTS OFF CACHE INTERNAL "")
set(AUTODIFF_BUILD_PYTHON OFF CACHE INTERNAL "")
set(AUTODIFF_BUILD_EXAMPLES OFF CACHE INTERNAL "")
set(AUTODIFF_BUILD_DOCS OFF CACHE INTERNAL "")

FetchContent_MakeAvailable( autodiff )

# Mark autodiff as system include directories so that no warnings trigger
get_target_property(_autodiff_includes autodiff INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(autodiff PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_autodiff_includes}")
