cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building Mille as part of the Acts project")

# Declare where to get Eigen from.
FetchContent_Declare(Mille SYSTEM ${ACTS_MILLE_SOURCE})

# Now set up its build.

# this prevents a target clash when Mille is built
# as child of another package
set(MILLE_AS_SUBDIR True)

FetchContent_MakeAvailable(Mille)
