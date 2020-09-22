# This file is part of the Acts project.
#
# Copyright (C) 2020 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# This module is used to find out whether the C++ compiler set up to build the
# project is able to compile SYCL code. If so, it provides helper variables
# to the project configuration to set up the build of SYCL
# libraries/executables.
#
# Variable(s) set up by the module:
#  - SYCL_FOUND: Boolean, set to TRUE when SYCL compilation is available.
#  - acts_target_setup_sycl( <target> [DEPENDENCY PUBLIC|PRIVATE|INTERFACE] ):
#       Helper function for setting up a target (library or executable) for
#       "SYCL compilation". Provides an optional argument for specifying
#       how/whether the target exposes this SYCL dependency.
#
# Technical/internal variable(s) set up by the module:
#  - SYCL_builtin_FOUND: Boolean, set to TRUE when the C++ compiler can itself
#                        understand/build SYCL code.
#  - SYCL_builtin_TARGETS: List of strings with the target platforms that the
#                          compiler can generate code for.
#  - SYCL_INCLUDE_DIR: Directory holding the (try)SYCL header(s).
#
# Variable(s) used to control the module's behaviour:
#  - SYCL_FLAGS: A list of compiler options that triggers a "SYCL build" from
#                the configured C++ compiler. Can be used to pass extra flags
#                to the SYCL compilation, if needed.
#  - SYCL_SEARCH_SUPPORT_LIBRARIES: A list of library (object file) names (more
#       like "name fragments") to include in linking SYCL binaries. By default
#       the code only looks for libsycl-glibc.o (or libsycl-msvc.o on Windows),
#       and libsycl-cmath.o. But one may need to use other object files (for
#       complex or 64-bit mathematics) coming with the Intel compiler as well.
#
# A typical use of the module would look like:
#
#  find_package( SYCL )
#  if( SYCL_FOUND )
#     atlas_add_library( SYCLAidedLibrary ... )
#     acts_target_setup_sycl( SYCLAidedLibrary DEPENDENCY PRIVATE )
#  endif()
#

# This module needs at least CMake 3.13.
cmake_minimum_required( VERSION 3.13 )

# We use CMake's built-in modules, used to look for SYCL capabilities.
include( CheckIncludeFileCXX )
include( CheckCXXSourceCompiles )
include( CMakeParseArguments )
set( CMAKE_REQUIRED_QUIET TRUE )

# Greet the user.
if( NOT SYCL_FIND_QUIETLY )
   message( STATUS "Checking if ${CMAKE_CXX_COMPILER} is SYCL capable..." )
endif()

# Set up the flag(s) needed for SYCL code compilation.
set( SYCL_FLAGS "-fsycl" CACHE STRING
   "Compiler flags to use for SYCL compilation" )
mark_as_advanced( SYCL_FLAGS )

# First check if the compiler is able to compile code using <CL/sycl.hpp> on
# its own, without any additional headers.
check_include_file_cxx( "CL/sycl.hpp" SYCL_builtin_FOUND "${SYCL_FLAGS}" )

# If that worked, we must be using a Clang version that understands sycl
# natively.
if( SYCL_builtin_FOUND )
   # Mark that SYCL is found.
   if( NOT SYCL_FIND_QUIETLY )
      message( STATUS
         "Checking if ${CMAKE_CXX_COMPILER} is SYCL capable... success" )
      message( STATUS "Checking for available SYCL target(s)..." )
   endif()
   set( SYCL_FOUND TRUE )
   # Figure out which SYCL target platforms are available.
   set( SYCL_POSSIBLE_TARGETS "spir64-unknown-unknown-sycldevice"
      CACHE STRING "List of targets to check the availability of" )
   mark_as_advanced( SYCL_POSSIBLE_TARGETS )
   set( SYCL_builtin_TARGETS )
   foreach( _target ${SYCL_POSSIBLE_TARGETS} )
      set( CMAKE_REQUIRED_FLAGS "-fsycl -fsycl-targets=${_target}" )
      check_cxx_source_compiles( "
         #include <CL/sycl.hpp>
         int main() {
            cl::sycl::platform::get_platforms();
            return 0;
         }
         " _syclTarget${_target}Found )
      if( _syclTarget${_target}Found )
         if( NOT SYCL_FIND_QUIETLY )
            message( STATUS "  - Found target: ${_target}" )
         endif()
         list( APPEND SYCL_builtin_TARGETS ${_target} )
      endif()
      unset( _syclTarget${_target}Found )
   endforeach()
   if( NOT SYCL_FIND_QUIETLY )
      message( STATUS "Checking for available SYCL target(s)... done" )
   endif()
   if( NOT "${SYCL_builtin_TARGETS}" STREQUAL "" )
      string( REPLACE ";" "," _targets "${SYCL_builtin_TARGETS}" )
      list( APPEND SYCL_FLAGS "-fsycl-targets=${_targets}" )
      unset( _targets )
   endif()
   # Look for object files holding SYCL device code, which would be needed for
   # the final binaries.
   set( SYCL_SEARCH_SUPPORT_LIBRARIES "" CACHE STRING
      "List of support libraries / object files to look for and link" )
   mark_as_advanced( SYCL_SEARCH_SUPPORT_LIBRARIES )
   get_filename_component( _compilerDir "${CMAKE_CXX_COMPILER}" DIRECTORY )
   set( SYCL_SUPPORT_LIBRARIES )
   foreach( _supportLib ${SYCL_SEARCH_SUPPORT_LIBRARIES} )
      find_file( SYCL_${_supportLib}_OBJECT_FILE
         NAMES "libsycl-${_supportLib}${CMAKE_CXX_OUTPUT_EXTENSION}"
         PATHS "${_compilerDir}"
         PATH_SUFFIXES "../lib" "../lib64" )
      if( SYCL_${_supportLib}_OBJECT_FILE )
         list( APPEND SYCL_SUPPORT_LIBRARIES
            "${SYCL_${_supportLib}_OBJECT_FILE}" )
      endif()
   endforeach()
   unset( _compilerDir )
   # Set up the acts_target_setup_sycl function.
   if( NOT COMMAND acts_target_setup_sycl )
      function( acts_target_setup_sycl targetName )
         cmake_parse_arguments( ARG "" "DEPENDENCY" "" ${ARGN} )
         if( NOT ARG_DEPENDENCY )
            set( ARG_DEPENDENCY "PRIVATE" )
         endif()
         target_compile_options( ${targetName} ${ARG_DEPENDENCY} ${SYCL_FLAGS} )
         target_link_options( ${targetName} ${ARG_DEPENDENCY} ${SYCL_FLAGS} )
         target_sources( ${targetName} ${ARG_DEPENDENCY}
            ${SYCL_SUPPORT_LIBRARIES} )
      endfunction( acts_target_setup_sycl )
   endif()
else()
   if( NOT SYCL_FIND_QUIETLY )
      message( STATUS
         "Checking if ${CMAKE_CXX_COMPILER} is SYCL capable... failure" )
      message( STATUS "Looking for the (tri)SYCL header(s)..." )
   endif()
   # If the compiler is not providing SYCL capabilities itself, check if we
   # have the (tri)SYCL headers available somewhere.
   find_path( SYCL_INCLUDE_DIR
      NAMES "CL/sycl.hpp" "SYCL/sycl.hpp" "triSYCL/sycl.hpp"
      PATH_SUFFIXES include ${CMAKE_INSTALL_INCLUDEDIR} include/triSYCL
      ${CMAKE_INSTALL_INCLUDEDIR}/triSYCL
      DOC "Location of the (tri)SYCL header(s)" )
   # Look for Boost as well, as these headers need it.
   if( SYCL_FIND_QUIETLY )
      find_package( Boost QUIET )
   else()
      find_package( Boost )
   endif()
   # Check if we found everything.
   if( SYCL_INCLUDE_DIR AND Boost_FOUND )
      # Mark that SYCL is found.
      if( NOT SYCL_FIND_QUIETLY )
         message( STATUS "Found (tri)SYCL headers: ${SYCL_INCLUDE_DIR}" )
      endif()
      set( SYCL_FOUND TRUE )
      # Set up the acts_target_setup_sycl function.
      if( NOT COMMAND acts_target_setup_sycl )
         function( acts_target_setup_sycl targetName )
            cmake_parse_arguments( ARG "" "DEPENDENCY" "" ${ARGN} )
            if( NOT ARG_DEPENDENCY )
               set( ARG_DEPENDENCY "PRIVATE" )
            endif()
            target_include_directories( ${targetName} SYSTEM ${ARG_DEPENDENCY}
               ${Boost_INCLUDE_DIRS} ${SYCL_INCLUDE_DIR} )
            target_link_libraries( ${targetName} ${ARG_DEPENDENCY}
               ${Boost_LIBRARIES} )
         endfunction( acts_target_setup_sycl )
      endif()
   else()
      # We did not find a viable SYCL version.
      if( NOT SYCL_FIND_QUIETLY )
         message( STATUS "Looking for the (tri)SYCL header(s)... failure" )
      endif()
      set( SYCL_FOUND FALSE )
   endif()
endif()
