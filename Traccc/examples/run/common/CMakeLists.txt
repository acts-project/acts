# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include( traccc-compiler-options-cpp )

# Create the common library.
traccc_add_library(traccc_examples_common examples_common TYPE SHARED
   "include/traccc/examples/make_magnetic_field.hpp"
   "src/make_magnetic_field.cpp"
   "include/traccc/examples/print_fitted_tracks_statistics.hpp"
   "src/print_fitted_tracks_statistics.cpp"
   "include/traccc/examples/throughput_mt.hpp"
   "include/traccc/examples/throughput_st.hpp"
   "include/traccc/examples/impl/throughput_mt.ipp"
   "include/traccc/examples/impl/throughput_st.ipp")
target_link_libraries(traccc_examples_common
   PUBLIC traccc::core traccc::options
   PRIVATE traccc::io)
