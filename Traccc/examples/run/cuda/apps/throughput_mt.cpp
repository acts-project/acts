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

int main(int argc, char* argv[]) {
  // We use a pinned CUDA host memory resource to allocate memory for the
  // cell inputs in order to speed up the memory copies.
  vecmem::cuda::host_memory_resource pinned_host_mr;

  // Execute the throughput test.
  return traccc::throughput_mt<traccc::cuda::full_chain_algorithm>(
      "Multi-threaded CUDA GPU throughput tests", argc, argv, &pinned_host_mr);
}
