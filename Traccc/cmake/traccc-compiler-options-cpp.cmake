# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Include the helper function(s).
include( traccc-functions )
include( CheckCXXCompilerFlag )

# Turn on the correct setting for the __cplusplus macro with MSVC.
if( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )
   traccc_add_flag( CMAKE_CXX_FLAGS "/Zc:__cplusplus" )
endif()

# Turn on a number of warnings for the "known compilers".
if( ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" ) OR
    ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) OR
    ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM" ) )

   # Basic flags for all build modes.
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wall" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wextra" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wshadow" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wunused-local-typedefs" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wpedantic" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wold-style-cast" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Wzero-as-null-pointer-constant" )
   traccc_add_flag( CMAKE_CXX_FLAGS "-Woverloaded-virtual" )
   if(PROJECT_IS_TOP_LEVEL)
     traccc_add_flag( CMAKE_CXX_FLAGS "-Wconversion" )
   endif()

   # Fail on warnings, if asked for that behaviour.
   if( TRACCC_FAIL_ON_WARNINGS )
      traccc_add_flag( CMAKE_CXX_FLAGS "-Werror" )
   endif()

elseif( "${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC" )

   # Basic flags for all build modes.
   string( REGEX REPLACE "/W[0-9]" "" CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS}" )
   traccc_add_flag( CMAKE_CXX_FLAGS "/W4" )

   # Fail on warnings, if asked for that behaviour.
   if( TRACCC_FAIL_ON_WARNINGS )
      traccc_add_flag( CMAKE_CXX_FLAGS "/WX" )
   endif()

endif()

# Set architecture level to v2 - enable generation of code using SSE4.2.
# But only use the flag if it's available from the compiler. And the user
# didn't specify some other -march flag yet.
if( "${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64" )
   set( X86_FLAG "-march=x86-64-v2" )
   check_cxx_compiler_flag( "${X86_FLAG}" TRACCC_SSE42_SUPPORTED )
   if( TRACCC_SSE42_SUPPORTED AND
       ( NOT "${CMAKE_CXX_FLAGS}" MATCHES "-march=" ) )
      traccc_add_flag( CMAKE_CXX_FLAGS "${X86_FLAG}" )
   endif()
   unset( X86_FLAG )
endif()
