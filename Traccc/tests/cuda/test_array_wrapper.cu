/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <traccc/utils/array_wrapper.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>

template <template <typename> typename F>
struct vector {
    using tuple_t = std::tuple<uint32_t, uint32_t, uint32_t>;
    F<uint32_t> x, y, z;
};

template <template <typename...> typename Layout>
__global__ void testArrayWrapperKernel(
    const typename traccc::array_wrapper<Layout, vector>::handle h,
    uint32_t *total) {
    __shared__ uint32_t block_total;

    if (threadIdx.x == 0) {
        block_total = 0;
    }

    __syncthreads();

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    uint32_t warp_total = h[tid].x;

    for (int i = 16; i >= 1; i /= 2) {
        warp_total += __shfl_xor_sync(0xffffffff, warp_total, i, 32);
    }

    if (threadIdx.x % 32 == 0) {
        atomicAdd(&block_total, warp_total);
    }

    __syncthreads();

    if (threadIdx.x == 0) {
        atomicAdd(total, block_total);
    }
}

template <template <typename...> typename Layout>
__global__ void fillWrapperKernel(
    typename traccc::array_wrapper<Layout, vector>::handle h) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < h.size()) {
        h[tid].x = 1u;
        h[tid].y = 2u;
        h[tid].z = 3u;
    }
}

TEST(CUDAArrayWrapper, SoALayout) {
    vecmem::cuda::device_memory_resource mr;

    uint32_t n = 1024u * 1024u;

    traccc::array_wrapper<traccc::soa, vector>::owner o(mr, n);

    fillWrapperKernel<traccc::soa><<<n / 1024u, 1024u>>>(
        typename traccc::array_wrapper<traccc::soa, vector>::handle(o));

    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    uint32_t host_result;
    uint32_t *dev_result = nullptr;

    ASSERT_EQ(cudaMalloc(&dev_result, sizeof(uint32_t)), cudaSuccess);
    ASSERT_EQ(cudaMemset(&dev_result, sizeof(uint32_t), 0), cudaSuccess);

    testArrayWrapperKernel<traccc::soa><<<n / 1024u, 1024u>>>(
        typename traccc::array_wrapper<traccc::soa, vector>::handle(o),
        dev_result);

    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    ASSERT_EQ(cudaMemcpy(&host_result, dev_result, sizeof(uint32_t),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);

    ASSERT_EQ(host_result, n);

    ASSERT_EQ(cudaFree(dev_result), cudaSuccess);
}

TEST(CUDAArrayWrapper, AoSLayout) {
    vecmem::cuda::device_memory_resource mr;

    uint32_t n = 1024u * 1024u;

    traccc::array_wrapper<traccc::aos, vector>::owner o(mr, n);

    fillWrapperKernel<traccc::aos><<<n / 1024u, 1024u>>>(
        typename traccc::array_wrapper<traccc::aos, vector>::handle(o));

    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    uint32_t host_result;
    uint32_t *dev_result = nullptr;

    ASSERT_EQ(cudaMalloc(&dev_result, sizeof(uint32_t)), cudaSuccess);
    ASSERT_EQ(cudaMemset(&dev_result, sizeof(uint32_t), 0), cudaSuccess);

    testArrayWrapperKernel<traccc::aos><<<n / 1024u, 1024u>>>(
        typename traccc::array_wrapper<traccc::aos, vector>::handle(o),
        dev_result);

    ASSERT_EQ(cudaPeekAtLastError(), cudaSuccess);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    ASSERT_EQ(cudaMemcpy(&host_result, dev_result, sizeof(uint32_t),
                         cudaMemcpyDeviceToHost),
              cudaSuccess);

    ASSERT_EQ(host_result, n);

    ASSERT_EQ(cudaFree(dev_result), cudaSuccess);
}
