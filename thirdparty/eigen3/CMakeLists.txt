cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building Eigen3 as part of the Acts project")

# Make sure Eigen exports an Eigen3::Eigen target that can be installed as
# part of Acts' own export set.
set(EIGEN_BUILD_CMAKE_PACKAGE ON CACHE BOOL "" FORCE)
# Declare where to get Eigen from.
FetchContent_Declare(Eigen SYSTEM ${ACTS_EIGEN3_SOURCE})
FetchContent_MakeAvailable(Eigen)
