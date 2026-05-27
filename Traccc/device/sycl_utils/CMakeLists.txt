# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# We need the SYCL language in this directory.
enable_language( SYCL )

# Project include(s).
include( traccc-compiler-options-cpp )
include( traccc-compiler-options-sycl )

# Set up the build of the traccc::sycl_utils library.
traccc_add_library( traccc_sycl_utils sycl_utils TYPE SHARED
  "include/traccc/sycl/utils/queue_wrapper.hpp"
  "src/queue_wrapper.cpp"
  "include/traccc/sycl/utils/make_magnetic_field.hpp"
  "src/make_magnetic_field.sycl" )
target_link_libraries( traccc_sycl_utils
  PUBLIC traccc::core
  PRIVATE covfie::core covfie::sycl )
