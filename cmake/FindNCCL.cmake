# Find the NVIDIA NCCL include directory and library.
#
# This module defines the `NCCL` imported target that encodes all
# necessary information in its target properties.

find_library(
  NCCL_LIBRARY
  NAMES nccl
  HINTS ${NCCL_DIR}
  PATH_SUFFIXES lib lib32 lib64
  DOC "The NCCL library")
  
if(NOT NCCL_LIBRARY)
  message(FATAL_ERROR "NCCL library not found")
endif()

find_path(
  NCCL_INCLUDE_DIR
  NAMES nccl.h nccl_net.h
  PATH_SUFFIXES include
  HITS ${NCCL_DIR}
  DOC "The NCCL include directory")
  
if(NOT OnnxRuntime_INCLUDE_DIR)
  message(FATAL_ERROR "NCCL includes not found")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  NCCL
  REQUIRED_VARS NCCL_LIBRARY NCCL_INCLUDE_DIR)

add_library(NCCL::NCCL SHARED IMPORTED)
set_property(TARGET NCCL::NCCL PROPERTY IMPORTED_LOCATION ${NCCL_LIBRARY})
set_property(TARGET NCCL::NCCL PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${NCCL_INCLUDE_DIR})

mark_as_advanced(NCCL_FOUND NCCL_INCLUDE_DIR NCCL_LIBRARY)
