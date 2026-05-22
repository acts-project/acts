cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building Millepede-II as part of the Acts project")

# Declare where to get Millepede-II from.
FetchContent_Declare(
    millepedeII
    SYSTEM
    ${ACTS_MP2_SOURCE}
    PATCH_COMMAND git apply ${CMAKE_CURRENT_SOURCE_DIR}/mp2.patch
)

# Now set up its build.
FetchContent_MakeAvailable(millepedeII)
