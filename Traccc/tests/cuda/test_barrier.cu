/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>

#include "../../device/cuda/src/utils/barrier.hpp"

__global__ void testBarrierAnd(bool* out) {
    traccc::cuda::barrier bar;

    bool v;

    v = bar.blockAnd(false);
    if (threadIdx.x == 0) {
        out[0] = v;
    }

    v = bar.blockAnd(true);
    if (threadIdx.x == 0) {
        out[1] = v;
    }

    v = bar.blockAnd(threadIdx.x % 2 == 0);
    if (threadIdx.x == 0) {
        out[2] = v;
    }

    v = bar.blockAnd(threadIdx.x < 32);
    if (threadIdx.x == 0) {
        out[3] = v;
    }
}

TEST(CUDABarrier, BarrierAnd) {
    vecmem::cuda::managed_memory_resource mr;
    constexpr std::size_t n_bools = 4;

    vecmem::unique_alloc_ptr<bool[]> out =
        vecmem::make_unique_alloc<bool[]>(mr, n_bools);

    testBarrierAnd<<<1, 1024>>>(out.get());

    ASSERT_EQ(cudaGetLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    EXPECT_FALSE(out.get()[0]);
    EXPECT_TRUE(out.get()[1]);
    EXPECT_FALSE(out.get()[2]);
    EXPECT_FALSE(out.get()[3]);
}

__global__ void testBarrierOr(bool* out) {
    traccc::cuda::barrier bar;

    bool v;

    v = bar.blockOr(false);
    if (threadIdx.x == 0) {
        out[0] = v;
    }

    v = bar.blockOr(true);
    if (threadIdx.x == 0) {
        out[1] = v;
    }

    v = bar.blockOr(threadIdx.x % 2 == 0);
    if (threadIdx.x == 0) {
        out[2] = v;
    }

    v = bar.blockOr(threadIdx.x < 32);
    if (threadIdx.x == 0) {
        out[3] = v;
    }
}

TEST(CUDABarrier, BarrierOr) {
    vecmem::cuda::managed_memory_resource mr;
    constexpr std::size_t n_bools = 4;

    vecmem::unique_alloc_ptr<bool[]> out =
        vecmem::make_unique_alloc<bool[]>(mr, n_bools);

    testBarrierOr<<<1, 1024>>>(out.get());

    ASSERT_EQ(cudaGetLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    EXPECT_FALSE(out.get()[0]);
    EXPECT_TRUE(out.get()[1]);
    EXPECT_TRUE(out.get()[2]);
    EXPECT_TRUE(out.get()[3]);
}

__global__ void testBarrierCount(int* out) {
    traccc::cuda::barrier bar;

    int v;

    v = bar.blockCount(false);
    if (threadIdx.x == 0) {
        out[0] = v;
    }

    v = bar.blockCount(true);
    if (threadIdx.x == 0) {
        out[1] = v;
    }

    v = bar.blockCount(threadIdx.x % 2 == 0);
    if (threadIdx.x == 0) {
        out[2] = v;
    }

    v = bar.blockCount(threadIdx.x < 32);
    if (threadIdx.x == 0) {
        out[3] = v;
    }
}

TEST(CUDABarrier, BarrierCount) {
    vecmem::cuda::managed_memory_resource mr;
    constexpr std::size_t n_ints = 4;

    vecmem::unique_alloc_ptr<int[]> out =
        vecmem::make_unique_alloc<int[]>(mr, n_ints);

    testBarrierCount<<<1, 1024>>>(out.get());

    ASSERT_EQ(cudaGetLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    EXPECT_EQ(out.get()[0], 0);
    EXPECT_EQ(out.get()[1], 1024);
    EXPECT_EQ(out.get()[2], 512);
    EXPECT_EQ(out.get()[3], 32);
}
