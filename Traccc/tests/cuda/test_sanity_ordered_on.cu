/*
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// vecmem includes
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

// traccc includes
#include <traccc/definitions/qualifiers.hpp>

#include "../../device/cuda/src/sanity/ordered_on.cuh"

// GTest include(s).
#include <gtest/gtest.h>

struct int_lt_relation {
    TRACCC_HOST_DEVICE
    bool operator()(const int& a, const int& b) const { return a < b; }
};

struct int_leq_relation {
    TRACCC_HOST_DEVICE
    bool operator()(const int& a, const int& b) const { return a <= b; }
};

class CUDASanityOrderedOn : public testing::Test {
    protected:
    vecmem::cuda::device_memory_resource mr;
    traccc::cuda::stream stream;
    vecmem::cuda::async_copy copy{stream.cudaStream()};
};

TEST_F(CUDASanityOrderedOn, TrueConsecutiveNoRepeatsLeq) {
    std::vector<int> host_vector;

    for (int i = 0; i < 500000; ++i) {
        host_vector.push_back(i);
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_TRUE(traccc::cuda::is_ordered_on<vecmem::device_vector<const int>>(
        int_leq_relation(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityOrderedOn, TrueConsecutiveNoRepeatsLt) {
    std::vector<int> host_vector;

    for (int i = 0; i < 500000; ++i) {
        host_vector.push_back(i);
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_TRUE(traccc::cuda::is_ordered_on<vecmem::device_vector<const int>>(
        int_lt_relation(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityOrderedOn, TrueConsecutiveRepeatsLeq) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_TRUE(traccc::cuda::is_ordered_on<vecmem::device_vector<const int>>(
        int_leq_relation(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityOrderedOn, FalseConsecutiveRepeatLt) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_FALSE(traccc::cuda::is_ordered_on<vecmem::device_vector<const int>>(
        int_lt_relation(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityOrderedOn, TrueConsecutivePathologicalFirstLeq) {
    std::vector<int> host_vector;

    host_vector.push_back(4000);

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_FALSE(traccc::cuda::is_ordered_on<vecmem::device_vector<const int>>(
        int_leq_relation(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityOrderedOn, TrueConsecutivePathologicalLastLeq) {
    std::vector<int> host_vector;

    host_vector.push_back(2000);

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_FALSE(traccc::cuda::is_ordered_on<vecmem::device_vector<const int>>(
        int_leq_relation(), mr, copy, stream, device_view));
}
