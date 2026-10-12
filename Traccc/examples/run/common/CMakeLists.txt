# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Project include(s).
include(traccc-compiler-options-cpp)

# Create the common library.
traccc_add_library(traccc_examples_common examples_common TYPE SHARED
   "include/traccc/examples/make_magnetic_field.hpp"
   "src/make_magnetic_field.cpp"
   "include/traccc/examples/await_strategy.hpp"
   "include/traccc/examples/print_fitted_tracks_statistics.hpp"
   "src/print_fitted_tracks_statistics.cpp"
   "include/traccc/examples/throughput_mt.hpp"
   "include/traccc/examples/throughput_st.hpp"
   "include/traccc/examples/impl/throughput_mt.ipp"
   "include/traccc/examples/impl/throughput_st.ipp"
)
target_link_libraries(
    traccc_examples_common
    PUBLIC traccc::core traccc::options
    PRIVATE traccc::io
)
