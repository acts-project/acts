/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/throughput_st.hpp"

#include "traccc/examples/hip/full_chain_algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/hip/host_memory_resource.hpp>

int main(int argc, char* argv[]) {
  // We use a pinned HIP host memory resource to allocate memory for the
  // cell inputs in order to speed up the memory copies.
  vecmem::hip::host_memory_resource pinned_host_mr;

  // Execute the throughput test.
  return traccc::throughput_st<traccc::hip::full_chain_algorithm>(
      "Single-threaded HIP GPU throughput tests", argc, argv, &pinned_host_mr);
}
