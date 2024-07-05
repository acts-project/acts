# Acts compiler flags
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type configuration" FORCE)
  message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif()

set(cxx_flags "-Wall -Wextra -Wpedantic -Wshadow -Wno-unused-local-typedefs")

# This adds some useful conversion checks like float-to-bool, float-to-int, etc.
# However, at the moment this is only added to clang builds, since GCC's -Wfloat-conversion
# is much more aggressive and also triggers on e.g., double-to-float
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(cxx_flags "${cxx_flags} -Wfloat-conversion")
endif()


set(ACTS_CXX_STANDARD 17)
set(ACTS_CXX_STANDARD_FEATURE cxx_std_17)
if(DEFINED CMAKE_CXX_STANDARD)
  if(${CMAKE_CXX_STANDARD} GREATER_EQUAL 17)
    set(ACTS_CXX_STANDARD ${CMAKE_CXX_STANDARD})
    set(ACTS_CXX_STANDARD_FEATURE "cxx_std_${CMAKE_CXX_STANDARD}")
  else()
    message(ERROR "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}, but ACTS requires C++ >=17")
  endif()
endif()


if(ACTS_ENABLE_CPU_PROFILING OR ACTS_ENABLE_MEMORY_PROFILING)
  message(STATUS "Added debug symbol compile flag")
  set(cxx_flags "${cxx_flags} ${CMAKE_CXX_FLAGS_DEBUG_INIT}")
endif()

# assign to global CXX flags
set(CMAKE_CXX_FLAGS "${cxx_flags} ${CMAKE_CXX_FLAGS}")
message(STATUS "Using compiler flags: ${CMAKE_CXX_FLAGS}")

# silence warning about missing RPATH on Mac OSX
set(CMAKE_MACOSX_RPATH 1)

# bake where we found external dependencies, if they were not in the default library directories
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# set relative library path for ACTS libraries
set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
