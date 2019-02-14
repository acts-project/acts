# set Acts compiler flags
set (ACTS_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wshadow -Wunused-local-typedefs")
set (ACTS_CXX_FLAGS_DEBUG "--coverage")
set (ACTS_CXX_FLAGS_MINSIZEREL "")
set (ACTS_CXX_FLAGS_RELEASE "")
set (ACTS_CXX_FLAGS_RELWITHDEBINFO "")

# set Acts linker flags
set (ACTS_EXE_LINKER_FLAGS_DEBUG "--coverage")
set (ACTS_SHARED_LINKER_FLAGS_DEBUG "--coverage ")

# assign to global CXX flags
set (CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${ACTS_CXX_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} ${ACTS_CXX_FLAGS_DEBUG}")
set (CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${ACTS_CXX_FLAGS_MINSIZEREL}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ACTS_CXX_FLAGS_RELEASE}")
set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${ACTS_CXX_FLAGS_RELWITHDEBINFO}")

# assign to global linker flags
set (CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${ACTS_EXE_LINKER_FLAGS_DEBUG}")
set (CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} ${ACTS_SHARED_LINKER_FLAGS_DEBUG}")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard version to use for the build")
if(${CMAKE_CXX_STANDARD} LESS 17)
  message(FATAL_ERROR "Acts cannot be build with a standard version\
  below C++17, C++${CMAKE_CXX_STANDARD} was requested")
endif()
message(STATUS "Building with standard version: C++${CMAKE_CXX_STANDARD}")

set(CMAKE_CXX_STANDARD_REQUIRED TRUE CACHE BOOL "Enforce C++ standard version.")
set(CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL "Allow/Disallow compiler extensions")

# silence warning about missing RPATH on Mac OSX
set (CMAKE_MACOSX_RPATH 1)
