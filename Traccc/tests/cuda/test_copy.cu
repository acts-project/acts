/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Projection include(s).
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/device/container_h2d_copy_alg.hpp"
#include "traccc/edm/container.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/fill.h>

// Google test include(s).
#include <gtest/gtest.h>

namespace {

using test_container_types = traccc::container_types<int, int>;

// Memory resources used by the application.
vecmem::host_memory_resource host_mr;
vecmem::cuda::device_memory_resource device_mr;
traccc::memory_resource mr{device_mr, &host_mr};

// Copy objects
vecmem::cuda::copy copy;
traccc::device::container_h2d_copy_alg<test_container_types> h2d{mr, copy};
traccc::device::container_d2h_copy_alg<test_container_types> d2h{mr, copy};

}  // namespace

TEST(CUDAContainerCopy, DeviceToHost) {

    // Create cuda container
    test_container_types::buffer device_buffer{{3, mr.main},
                                               {{1, 3, 2}, mr.main, mr.host}};

    // Get device vectors and fill them in the device
    vecmem::device_vector<int> device_headers(device_buffer.headers);
    vecmem::device_vector<int> device_items_0(
        device_buffer.items.host_ptr()[0]);
    vecmem::device_vector<int> device_items_1(
        device_buffer.items.host_ptr()[1]);
    vecmem::device_vector<int> device_items_2(
        device_buffer.items.host_ptr()[2]);

    thrust::fill(thrust::device, device_headers.begin(), device_headers.end(),
                 7);
    thrust::fill(thrust::device, device_items_0.begin(), device_items_0.end(),
                 4);
    thrust::fill(thrust::device, device_items_1.begin(), device_items_1.end(),
                 1);
    thrust::fill(thrust::device, device_items_2.begin(), device_items_2.end(),
                 2);

    // Device-to-Host Copy
    test_container_types::host host_copy = d2h(device_buffer);

    // Check the copied container
    const auto& host_headers = host_copy.get_headers();
    ASSERT_EQ(host_headers.size(), 3u);
    ASSERT_EQ(host_headers[0], 7);
    ASSERT_EQ(host_headers[1], 7);
    ASSERT_EQ(host_headers[2], 7);

    const auto& host_items = host_copy.get_items();
    ASSERT_EQ(host_items.size(), 3u);
    ASSERT_EQ(host_items[0].size(), 1u);
    ASSERT_EQ(host_items[0][0], 4);
    ASSERT_EQ(host_items[1].size(), 3u);
    ASSERT_EQ(host_items[1][0], 1);
    ASSERT_EQ(host_items[1][1], 1);
    ASSERT_EQ(host_items[1][2], 1);
    ASSERT_EQ(host_items[2].size(), 2u);
    ASSERT_EQ(host_items[2][0], 2);
    ASSERT_EQ(host_items[2][1], 2);
}

TEST(CUDAContainerCopy, HostToDeviceToHost) {

    // Create a host container
    test_container_types::host host_orig{&host_mr};

    // Fill the host container
    for (int i = 0; i < 3; ++i) {
        host_orig.push_back(
            i, test_container_types::host::vector_type<int>{i, i + 1});
    }

    // Copy the container to the device, and back.
    test_container_types::host host_copy =
        d2h(h2d(traccc::get_data(host_orig)));

    // Check the copied container
    ASSERT_EQ(host_copy.size(), 3u);
    const auto& host_headers = host_copy.get_headers();
    ASSERT_EQ(host_headers.size(), 3u);
    ASSERT_EQ(host_headers[0], 0);
    ASSERT_EQ(host_headers[1], 1);
    ASSERT_EQ(host_headers[2], 2);

    const auto& host_items = host_copy.get_items();
    ASSERT_EQ(host_items.size(), 3u);
    ASSERT_EQ(host_items[0][0], 0);
    ASSERT_EQ(host_items[0][1], 1);
    ASSERT_EQ(host_items[1][0], 1);
    ASSERT_EQ(host_items[1][1], 2);
    ASSERT_EQ(host_items[2][0], 2);
    ASSERT_EQ(host_items[2][1], 3);
}
