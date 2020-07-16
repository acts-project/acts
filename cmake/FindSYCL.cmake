# Greet the user.
if( NOT SYCL_FIND_QUIETLY )
   message( STATUS "Checking if ${CMAKE_CXX_COMPILER} is SYCL capable..." )
endif()

# First check if the compiler is able to compile code using <CL/sycl.hpp> on
# its own, without any additional headers.
check_include_file_cxx( "CL/sycl.hpp" SYCL_builtin_FOUND "-fsycl" )

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
   set( SYCL_FLAGS "-fsycl" )
   set( SYCL_builtin_TARGETS )
   foreach( _target "nvptx64-*-*-sycldevice" "spir64-*-*-sycldevice" )
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
   # Set up the atlas_setup_sycl_target function.
   if( NOT COMMAND atlas_setup_sycl_target )
      function( atlas_setup_sycl_target targetName )
         cmake_parse_arguments( ARG "" "DEPENDENCY" "" ${ARGN} )
         if( NOT ARG_DEPENDENCY )
            set( ARG_DEPENDENCY "PRIVATE" )
         endif()
         target_compile_options( ${targetName} ${ARG_DEPENDENCY} ${SYCL_FLAGS} )
         target_link_options( ${targetName} ${ARG_DEPENDENCY} ${SYCL_FLAGS} )
      endfunction( atlas_setup_sycl_target )
   endif()
else()
   # We did not find a viable SYCL version.
   if( NOT SYCL_FIND_QUIETLY )
   message( STATUS
      "Checking if ${CMAKE_CXX_COMPILER} is SYCL capable... failure" )
   endif()
   set( SYCL_FOUND FALSE )
endif()