# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include( traccc-functions )

# Warning flags for the AMD backend of the HIP compiler.
if( "${CMAKE_HIP_PLATFORM}" STREQUAL "amd" )
   traccc_add_flag( CMAKE_HIP_FLAGS "-Wall" )
   traccc_add_flag( CMAKE_HIP_FLAGS "-Wextra" )
   traccc_add_flag( CMAKE_HIP_FLAGS "-Wshadow" )
   traccc_add_flag( CMAKE_HIP_FLAGS "-Wunused-local-typedefs" )
   traccc_add_flag( CMAKE_HIP_FLAGS "-pedantic" )
endif()

# Specific flags for the NVIDIA backend of the HIP compiler.
if( "${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia" )
   traccc_add_flag( CMAKE_HIP_FLAGS "--expt-relaxed-constexpr" )
   traccc_add_flag( CMAKE_HIP_FLAGS "--use_fast_math" )
endif()

# Fail on warnings, if asked for that behaviour.
if( TRACCC_FAIL_ON_WARNINGS )
   if( "${CMAKE_HIP_PLATFORM}" STREQUAL "amd" )
      traccc_add_flag( CMAKE_HIP_FLAGS "-Werror" )
   elseif( "${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia" )
      traccc_add_flag( CMAKE_HIP_FLAGS "-Werror all-warnings" )
   endif()
endif()
