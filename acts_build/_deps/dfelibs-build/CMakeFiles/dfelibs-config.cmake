# DFE library, part of the ACTS project
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set up the helper functions/macros.

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was dfelibs-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

# Set up some simple variables for using the package.
set(dfelibs_VERSION "20200416")
set_and_check(dfelibs_INCLUDE_DIR
  "${PACKAGE_PREFIX_DIR}/include")
set_and_check(dfelibs_CMAKE_DIR "${PACKAGE_PREFIX_DIR}/lib/cmake/dfelibs-20200416")

# Print a standard information message about the package being found.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(dfelibs REQUIRED_VARS
   CMAKE_CURRENT_LIST_FILE
   VERSION_VAR dfelibs_VERSION)

# Include the file listing all the imported targets and options.
include("${dfelibs_CMAKE_DIR}/dfelibs-config-targets.cmake")
