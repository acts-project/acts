# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set the CPU build flags.
include(detray-compiler-options-cpp)

# Set up a test library, which the "new style" benchmarks and tests could use.
add_library(
    detray_benchmark_cpu
    INTERFACE
    "matrix_benchmark.hpp"
    "vector_benchmark.hpp"
    "transform_benchmark.hpp"
    "propagation_benchmark.hpp"
)

add_library(detray::benchmark_cpu ALIAS detray_benchmark_cpu)

target_link_libraries(detray_benchmark_cpu INTERFACE detray::benchmarks)
