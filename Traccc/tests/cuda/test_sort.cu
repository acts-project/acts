/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>

#include "../../cuda/src/utils/barrier.hpp"
#include "../../cuda/src/utils/thread_id.hpp"
#include "traccc/device/sort.hpp"

__global__ void testBlockSortKernel(std::uint32_t *keys, std::uint32_t n_keys) {
  traccc::cuda::details::thread_id1 thread_id;
  traccc::cuda::barrier barrier;
  traccc::device::blockOddEvenSort(thread_id, barrier, keys, n_keys,
                                   std::less<std::uint32_t>());
}

TEST(CUDASort, BlockOddEvenSort) {
  vecmem::cuda::managed_memory_resource mr;

  std::uint32_t n = 2803;
  vecmem::unique_alloc_ptr<std::uint32_t[]> arr =
      vecmem::make_unique_alloc<std::uint32_t[]>(mr, n);

  // As long as 13 and n_keys are coprime, this will generate a big,
  // non-sorted array containing every element.
  for (std::uint32_t i = 0; i < n; i++) {
    arr[i] = (13 * 500 * i) % n;
  }

  testBlockSortKernel<<<1, 1024u>>>(arr.get(), n);

  ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
  ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

  for (std::uint32_t i = 0; i < n; ++i) {
    ASSERT_EQ(arr[i], i);
  }
}
