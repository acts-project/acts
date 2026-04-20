cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building Millepede-II as part of the Acts project")

# Declare where to get Eigen from.
FetchContent_Declare(millepedeII SYSTEM ${ACTS_MP2_SOURCE})

# Now set up its build.
FetchContent_MakeAvailable(millepedeII)
