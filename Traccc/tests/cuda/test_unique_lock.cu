/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <mutex>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>

#include "../../device/cuda/src/utils/cuda_error_handling.hpp"
#include "traccc/device/mutex.hpp"
#include "traccc/device/unique_lock.hpp"

__global__ void unique_lock_add_kernel_try_lock(uint32_t *out,
                                                uint32_t *_lock) {
    traccc::device::mutex m(*_lock);

    if (threadIdx.x == 0) {
        traccc::device::unique_lock lock(m, std::try_to_lock);

        if (!lock.owns_lock()) {
            lock.lock();
        }

        uint32_t tmp = *out;
        tmp += 1;
        *out = tmp;
    }
}

__global__ void unique_lock_add_kernel_defer_lock(uint32_t *out,
                                                  uint32_t *_lock) {
    traccc::device::mutex m(*_lock);
    traccc::device::unique_lock lock(m, std::defer_lock);

    if (threadIdx.x == 0) {
        lock.lock();

        uint32_t tmp = *out;
        tmp += 1;
        *out = tmp;
    }
}

__global__ void unique_lock_add_kernel_adopt_lock(uint32_t *out,
                                                  uint32_t *_lock) {
    traccc::device::mutex m(*_lock);

    if (threadIdx.x == 0) {
        m.lock();
        traccc::device::unique_lock lock(m, std::adopt_lock);

        uint32_t tmp = *out;
        tmp += 1;
        *out = tmp;
    }
}

TEST(CUDAUniqueLock, MassAdditionKernelTryLock) {
    vecmem::cuda::managed_memory_resource mr;

    vecmem::unique_alloc_ptr<uint32_t> out =
        vecmem::make_unique_alloc<uint32_t>(mr);
    vecmem::unique_alloc_ptr<uint32_t> lock =
        vecmem::make_unique_alloc<uint32_t>(mr);

    TRACCC_CUDA_ERROR_CHECK(cudaMemset(lock.get(), 0, sizeof(uint32_t)));
    TRACCC_CUDA_ERROR_CHECK(cudaMemset(out.get(), 0, sizeof(uint32_t)));

    uint32_t n_blocks = 262144;
    uint32_t n_threads = 32;

    unique_lock_add_kernel_try_lock<<<n_blocks, n_threads>>>(out.get(),
                                                             lock.get());

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    EXPECT_EQ(n_blocks, *out.get());
}

TEST(CUDAUniqueLock, MassAdditionKernelDeferLock) {
    vecmem::cuda::managed_memory_resource mr;

    vecmem::unique_alloc_ptr<uint32_t> out =
        vecmem::make_unique_alloc<uint32_t>(mr);
    vecmem::unique_alloc_ptr<uint32_t> lock =
        vecmem::make_unique_alloc<uint32_t>(mr);

    TRACCC_CUDA_ERROR_CHECK(cudaMemset(lock.get(), 0, sizeof(uint32_t)));
    TRACCC_CUDA_ERROR_CHECK(cudaMemset(out.get(), 0, sizeof(uint32_t)));

    uint32_t n_blocks = 262144;
    uint32_t n_threads = 32;

    unique_lock_add_kernel_defer_lock<<<n_blocks, n_threads>>>(out.get(),
                                                               lock.get());

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    EXPECT_EQ(n_blocks, *out.get());
}

TEST(CUDAUniqueLock, MassAdditionKernelAdoptLock) {
    vecmem::cuda::managed_memory_resource mr;

    vecmem::unique_alloc_ptr<uint32_t> out =
        vecmem::make_unique_alloc<uint32_t>(mr);
    vecmem::unique_alloc_ptr<uint32_t> lock =
        vecmem::make_unique_alloc<uint32_t>(mr);

    TRACCC_CUDA_ERROR_CHECK(cudaMemset(lock.get(), 0, sizeof(uint32_t)));
    TRACCC_CUDA_ERROR_CHECK(cudaMemset(out.get(), 0, sizeof(uint32_t)));

    uint32_t n_blocks = 262144;
    uint32_t n_threads = 32;

    unique_lock_add_kernel_adopt_lock<<<n_blocks, n_threads>>>(out.get(),
                                                               lock.get());

    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());

    EXPECT_EQ(n_blocks, *out.get());
}
