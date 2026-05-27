/*
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// vecmem includes
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>

// traccc includes
#include <traccc/definitions/qualifiers.hpp>
#include <traccc/sanity/contiguous_on.hpp>

// GTest include(s).
#include <gtest/gtest.h>

struct int_identity_projection {
    TRACCC_HOST_DEVICE
    int operator()(const int& v) const { return v; }
};

class CPUSanityContiguousOn : public testing::Test {
    protected:
    CPUSanityContiguousOn() {}
};

TEST_F(CPUSanityContiguousOn, TrueOrdered) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    ASSERT_TRUE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_TRUE(traccc::host::is_contiguous_on(int_identity_projection(),
                                               device_vector));
}

TEST_F(CPUSanityContiguousOn, TrueRandom) {
    std::vector<int> host_vector;

    for (int i : {603, 6432, 1, 3, 67, 2, 1111}) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    ASSERT_TRUE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_TRUE(traccc::host::is_contiguous_on(int_identity_projection(),
                                               device_vector));
}

TEST_F(CPUSanityContiguousOn, FalseOrdered) {
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

    ASSERT_FALSE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_FALSE(traccc::host::is_contiguous_on(int_identity_projection(),
                                                device_vector));
}

TEST_F(CPUSanityContiguousOn, FalseOrderedPathologicalFirst) {
    std::vector<int> host_vector;

    host_vector.push_back(4000);

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    ASSERT_FALSE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_FALSE(traccc::host::is_contiguous_on(int_identity_projection(),
                                                device_vector));
}

TEST_F(CPUSanityContiguousOn, TrueOrderedPathologicalFirst) {
    std::vector<int> host_vector;

    host_vector.push_back(6000);

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    ASSERT_TRUE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_TRUE(traccc::host::is_contiguous_on(int_identity_projection(),
                                               device_vector));
}

TEST_F(CPUSanityContiguousOn, FalseOrderedPathologicalLast) {
    std::vector<int> host_vector;

    for (int i = 0; i < 5000; ++i) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    host_vector.push_back(2);

    ASSERT_FALSE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_FALSE(traccc::host::is_contiguous_on(int_identity_projection(),
                                                device_vector));
}

TEST_F(CPUSanityContiguousOn, FalseRandom) {
    std::vector<int> host_vector;

    for (int i : {603, 6432, 1, 3, 67, 1, 1111}) {
        for (int j = 0; j < i; ++j) {
            host_vector.push_back(i);
        }
    }

    ASSERT_FALSE(
        traccc::host::is_contiguous_on(int_identity_projection(), host_vector));

    vecmem::device_vector<int> device_vector(vecmem::get_data(host_vector));

    ASSERT_FALSE(traccc::host::is_contiguous_on(int_identity_projection(),
                                                device_vector));
}
