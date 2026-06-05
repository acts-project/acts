/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "../../device/cuda/src/utils/sync.cuh"

__global__ void testWarpIndexedBallotSyncBasicKernel(uint32_t *vts,
                                                     uint32_t *vis) {
    auto [vt, vi] =
        traccc::cuda::warp_indexed_ballot_sync(threadIdx.x % 2 == 0);

    vts[threadIdx.x] = vt;
    vis[threadIdx.x] = vi;
}

__global__ void testWarpIndexedBallotSyncWithExitKernel(uint32_t *vts,
                                                        uint32_t *vis) {
    if (threadIdx.x < 16) {
        return;
    }

    auto [vt, vi] =
        traccc::cuda::warp_indexed_ballot_sync(threadIdx.x % 2 == 0);

    vts[threadIdx.x] = vt;
    vis[threadIdx.x] = vi;
}

TEST(CUDASync, WarpIndexedBallotSyncBasic) {
    uint32_t *dev_vt = nullptr, *dev_vi = nullptr;
    uint32_t host_vt[32], host_vi[32];

    ASSERT_EQ(cudaMalloc(&dev_vt, 32u * sizeof(uint32_t)), cudaSuccess);
    ASSERT_EQ(cudaMalloc(&dev_vi, 32u * sizeof(uint32_t)), cudaSuccess);
    ASSERT_NE(dev_vt, nullptr);
    ASSERT_NE(dev_vi, nullptr);

    testWarpIndexedBallotSyncBasicKernel<<<1, 32u>>>(dev_vt, dev_vi);

    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);

    ASSERT_EQ(cudaMemcpy(host_vt, dev_vt, 32u * sizeof(uint32_t),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);
    ASSERT_EQ(cudaMemcpy(host_vi, dev_vi, 32u * sizeof(uint32_t),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);

    for (uint32_t i = 0; i < 32u; ++i) {
        ASSERT_EQ(host_vt[i], 16u);
    }

    for (uint32_t i = 0; i < 16u; ++i) {
        ASSERT_EQ(host_vi[i * 2], i);
    }

    ASSERT_EQ(cudaFree(dev_vt), cudaSuccess);
    ASSERT_EQ(cudaFree(dev_vi), cudaSuccess);
}

TEST(CUDASync, WarpIndexedBallotSyncWithExit) {
    uint32_t *dev_vt = nullptr, *dev_vi = nullptr;
    uint32_t host_vt[32], host_vi[32];

    ASSERT_EQ(cudaMalloc(&dev_vt, 32u * sizeof(uint32_t)), cudaSuccess);
    ASSERT_EQ(cudaMalloc(&dev_vi, 32u * sizeof(uint32_t)), cudaSuccess);
    ASSERT_NE(dev_vt, nullptr);
    ASSERT_NE(dev_vi, nullptr);

    testWarpIndexedBallotSyncWithExitKernel<<<1, 32u>>>(dev_vt, dev_vi);

    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);

    ASSERT_EQ(cudaMemcpy(host_vt, dev_vt, 32u * sizeof(uint32_t),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);
    ASSERT_EQ(cudaMemcpy(host_vi, dev_vi, 32u * sizeof(uint32_t),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);

    for (uint32_t i = 16; i < 32u; ++i) {
        ASSERT_EQ(host_vt[i], 8u);
    }

    for (uint32_t i = 0; i < 8u; ++i) {
        ASSERT_EQ(host_vi[i * 2 + 16], i);
    }

    ASSERT_EQ(cudaFree(dev_vt), cudaSuccess);
    ASSERT_EQ(cudaFree(dev_vi), cudaSuccess);
}
