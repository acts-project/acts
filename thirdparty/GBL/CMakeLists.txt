cmake_minimum_required(VERSION 3.25)

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(
    STATUS
    "Building General Broken Lines (GBL) from Mille as part of the Acts project"
)

# Declare where to get GBL from.
FetchContent_Declare(GBL SYSTEM ${ACTS_GBL_SOURCE} SOURCE_SUBDIR cpp)

# Now set up its build.
FetchContent_MakeAvailable(GBL)
