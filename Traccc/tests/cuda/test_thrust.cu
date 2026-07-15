/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

// GTest include(s).
#include <gtest/gtest.h>

// This defines the local frame test suite

namespace {
vecmem::cuda::copy copy;
vecmem::host_memory_resource host_resource;
vecmem::cuda::device_memory_resource device_resource;

}  // namespace

TEST(CUDAThrust, sort) {

    vecmem::vector<unsigned int> host_vector{{3, 2, 1, 8, 4}, &host_resource};

    auto host_buffer = vecmem::get_data(host_vector);
    auto device_buffer = copy.to(vecmem::get_data(host_vector), device_resource,
                                 vecmem::copy::type::host_to_device);

    vecmem::device_vector<unsigned int> device_vector(device_buffer);

    thrust::sort(thrust::device, device_vector.begin(), device_vector.end());

    copy(device_buffer, host_buffer, vecmem::copy::type::device_to_host)
        ->wait();
    ;

    ASSERT_EQ(host_vector[0], 1);
    ASSERT_EQ(host_vector[1], 2);
    ASSERT_EQ(host_vector[2], 3);
    ASSERT_EQ(host_vector[3], 4);
    ASSERT_EQ(host_vector[4], 8);
}

TEST(CUDAThrust, scan) {

    vecmem::vector<unsigned int> host_vector{{3, 2, 1, 8, 4}, &host_resource};

    auto host_buffer = vecmem::get_data(host_vector);
    auto device_buffer = copy.to(vecmem::get_data(host_vector), device_resource,
                                 vecmem::copy::type::host_to_device);

    vecmem::device_vector<unsigned int> device_vector(device_buffer);

    thrust::inclusive_scan(thrust::device, device_vector.begin(),
                           device_vector.end(), device_vector.begin());

    copy(device_buffer, host_buffer, vecmem::copy::type::device_to_host)
        ->wait();

    ASSERT_EQ(host_vector[0], 3);
    ASSERT_EQ(host_vector[1], 5);
    ASSERT_EQ(host_vector[2], 6);
    ASSERT_EQ(host_vector[3], 14);
    ASSERT_EQ(host_vector[4], 18);
}

TEST(CUDAThrust, fill) {

    vecmem::vector<unsigned int> host_vector{{1, 1, 1, 1, 1, 1, 1},
                                             &host_resource};

    auto host_buffer = vecmem::get_data(host_vector);
    auto device_buffer = copy.to(vecmem::get_data(host_vector), device_resource,
                                 vecmem::copy::type::host_to_device);

    vecmem::device_vector<unsigned int> device_vector(device_buffer);

    thrust::fill(thrust::device, device_vector.begin(), device_vector.end(),
                 112);

    copy(device_buffer, host_buffer, vecmem::copy::type::device_to_host)
        ->wait();

    ASSERT_EQ(host_vector[0], 112);
    ASSERT_EQ(host_vector[1], 112);
    ASSERT_EQ(host_vector[2], 112);
    ASSERT_EQ(host_vector[3], 112);
    ASSERT_EQ(host_vector[4], 112);
    ASSERT_EQ(host_vector[5], 112);
    ASSERT_EQ(host_vector[6], 112);
}
