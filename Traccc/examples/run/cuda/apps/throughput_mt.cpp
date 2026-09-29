/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/throughput_mt.hpp"

#include "traccc/examples/cuda/full_chain_algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/host_memory_resource.hpp>

// System include(s).
#include <cstdlib>

int main(int argc, char* argv[]) {
  // We use a pinned CUDA host memory resource to allocate memory for the
  // cell inputs in order to speed up the memory copies.
  vecmem::cuda::host_memory_resource pinned_host_mr;

  // Using pinned host memory interferes with NSight Compute, so we try to
  // detect whether we are running inside ncu; if so, we use a non-pinned
  // MR.
  //
  // WARNING: This is based on undocumented environment variables that ncu
  // sets, so this may break (harmlessly) in the future.
  vecmem::memory_resource* const resource_ptr =
      (std::getenv("NV_COMPUTE_PROFILER_PERFWORKS_DIR") != nullptr)
          ? nullptr
          : &pinned_host_mr;

  // Execute the throughput test.
  return traccc::throughput_mt<traccc::cuda::full_chain_algorithm>(
      "Multi-threaded CUDA GPU throughput tests", argc, argv, resource_ptr);
}
