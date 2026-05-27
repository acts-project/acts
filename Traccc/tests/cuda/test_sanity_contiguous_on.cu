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
#include "../../device/cuda/src/sanity/contiguous_on.cuh"
#include "traccc/definitions/qualifiers.hpp"

// GTest include(s).
#include <gtest/gtest.h>

struct int_identity_projection {
    TRACCC_HOST_DEVICE
    int operator()(const int& v) const { return v; }
};

class CUDASanityContiguousOn : public testing::Test {
    protected:
    vecmem::cuda::device_memory_resource mr;
    traccc::cuda::stream stream;
    vecmem::cuda::async_copy copy{stream.cudaStream()};
};

TEST_F(CUDASanityContiguousOn, TrueOrdered) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_TRUE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityContiguousOn, TrueRandom) {
    std::vector<int> host_vector;

    for (int i : {603, 6432, 1, 3, 67, 2, 1111}) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_TRUE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityContiguousOn, FalseOrdered) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        if (i == 105) {
            host_vector.push_back(5);
        } else {
            for (int j = 0; j < i; ++j) {
                host_vector.push_back(i);
            }
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_FALSE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityContiguousOn, FalseOrderedPathologicalFirst) {
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

    ASSERT_FALSE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityContiguousOn, TrueOrderedPathologicalFirst) {
    std::vector<int> host_vector;

    host_vector.push_back(6000);

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_TRUE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityContiguousOn, FalseOrderedPathologicalLast) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    host_vector.push_back(2);

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_FALSE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}

TEST_F(CUDASanityContiguousOn, FalseRandom) {
    std::vector<int> host_vector;

    for (int i : {603, 6432, 1, 3, 67, 1, 1111}) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    auto device_data = copy.to(vecmem::get_data(host_vector), mr,
                               vecmem::copy::type::host_to_device);
    auto device_view = vecmem::get_data(device_data);

    ASSERT_FALSE(
        traccc::cuda::is_contiguous_on<vecmem::device_vector<const int>>(
            int_identity_projection(), mr, copy, stream, device_view));
}
