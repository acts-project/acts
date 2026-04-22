# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

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
