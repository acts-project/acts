/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"

// Traccc test include(s)
#include "tests/grid2_test.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>

using namespace traccc;

TEST(grids_cuda, grid2_replace_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis2::regular<> xaxis{4u, -1.f, 3.f, mng_mr};
    axis2::regular<> yaxis{6u, 0.f, 6.f, mng_mr};

    auto x_interval =
        (xaxis.max - xaxis.min) / static_cast<scalar>(xaxis.n_bins);
    auto y_interval =
        (yaxis.max - yaxis.min) / static_cast<scalar>(yaxis.n_bins);

    // declare host grid
    host_grid2_replace g2(std::move(xaxis), std::move(yaxis), mng_mr,
                          point3{0.f, 0.f, 0.f});

    // pre-check
    for (unsigned int i_x = 0u; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < yaxis.bins(); i_y++) {

            const auto& data = g2.bin(i_x, i_y);

            EXPECT_EQ(data, g2.populator().m_invalid);
        }
    }

    // get grid_data
    auto g2_data = get_data(g2, mng_mr);

    // fill the grids
    grid_replace_test(g2_data);

    // post-check
    for (unsigned int i_x = 0u; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < yaxis.bins(); i_y++) {
            auto bin_id = static_cast<scalar>(i_x + i_y * xaxis.bins());
            const auto& data = g2.bin(i_x, i_y);

            point3 tp({xaxis.min + bin_id * x_interval,
                       yaxis.min + bin_id * y_interval, 0.5f});

            EXPECT_EQ(data, tp);
        }
    }
}

TEST(grids_cuda, grid2_replace_populator_ci) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis2::circular<> caxis{4u, -2.f, 2.f, mng_mr};
    axis2::irregular<> iaxis{{1.f, 3.f, 9.f, 27.f, 81.f}, mng_mr};

    auto x_interval =
        (caxis.max - caxis.min) / static_cast<scalar>(caxis.n_bins);

    // declare host grid
    host_grid2_replace_ci g2(std::move(caxis), std::move(iaxis), mng_mr,
                             point3{0.f, 0.f, 0.f});

    // pre-check
    for (unsigned int i_x = 0u; i_x < caxis.bins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < iaxis.bins(); i_y++) {

            const auto& data = g2.bin(i_x, i_y);

            EXPECT_EQ(data, g2.populator().m_invalid);
        }
    }

    // get grid_data
    auto g2_data = get_data(g2, mng_mr);

    // fill the grids
    grid_replace_ci_test(g2_data);

    // post-check
    for (unsigned int i_x = 0u; i_x < caxis.bins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < iaxis.bins(); i_y++) {
            auto y_interval = iaxis.boundaries[i_y + 1] - iaxis.boundaries[i_y];
            auto bin_id = static_cast<scalar>(i_x + i_y * caxis.bins());

            const auto& data = g2.bin(i_x, i_y);

            point3 tp({caxis.min + bin_id * x_interval,
                       iaxis.min + bin_id * y_interval, 0.5f});

            EXPECT_EQ(data, tp);
        }
    }
}

// Equality operator in complete populator does not work correctly in CUDA
// (!constexpr)
TEST(grids_cuda, grid2_complete_populator) {
    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    // axis
    axis2::regular<> xaxis{7u, -1.f, 6.f, mng_mr};
    axis2::regular<> yaxis{3u, 0.f, 3.f, mng_mr};

    // declare grid
    host_grid2_complete g2(std::move(xaxis), std::move(yaxis), mng_mr,
                           point3{0.f, 0.f, 0.f});

    // pre-check
    for (unsigned int i_x = 0u; i_x < xaxis.bins(); i_x++) {
        for (unsigned int i_y = 0u; i_y < yaxis.bins(); i_y++) {

            const auto& data = g2.bin(i_x, i_y);

            for (auto pt : data) {
                EXPECT_EQ(pt, g2.populator().m_invalid);
            }
        }
    }

    // get grid_data
    auto g2_data = get_data(g2, mng_mr);

    // fill the grid
    grid_complete_test(g2_data);

    auto x_interval =
        (xaxis.max - xaxis.min) / static_cast<scalar>(xaxis.n_bins);
    auto y_interval =
        (yaxis.max - yaxis.min) / static_cast<scalar>(yaxis.n_bins);

    // post-check
    for (unsigned int i_y = 0u; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0u; i_x < xaxis.bins(); i_x++) {

            const auto& data = g2.bin(i_x, i_y);

            for (unsigned int i_p = 0u; i_p < data.size(); i_p++) {
                auto& pt = data[i_p];

                auto bin_id = i_x + i_y * xaxis.bins();
                auto gid = static_cast<scalar>(i_p + bin_id * data.size());

                point3 tp({xaxis.min + gid * x_interval,
                           yaxis.min + gid * y_interval, 0.5f});
                EXPECT_EQ(pt, tp);
            }
        }
    }
}

TEST(grids_cuda, grid2_attach_populator) {

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    axis2::circular<> xaxis{65u, -constant<scalar>::pi, constant<scalar>::pi,
                            mng_mr};
    axis2::regular<> yaxis{2u, 0.f, 6.f, mng_mr};

    auto x_interval =
        (xaxis.max - xaxis.min) / static_cast<scalar>(xaxis.n_bins);
    auto y_interval =
        (yaxis.max - yaxis.min) / static_cast<scalar>(yaxis.n_bins);

    host_grid2_attach g2(xaxis, yaxis, mng_mr, point3{0.f, 0.f, 0.f});

    for (unsigned int i_y = 0u; i_y < yaxis.bins(); i_y++) {
        for (unsigned int i_x = 0u; i_x < xaxis.bins(); i_x++) {

            for (unsigned int i_p = 0u; i_p < 100u; i_p++) {

                auto bin_id = i_x + i_y * xaxis.bins();
                auto gid = static_cast<scalar>(i_p + bin_id * 100u);

                point3 tp({xaxis.min + gid * x_interval,
                           yaxis.min + gid * y_interval, 0.5f});
                g2.populate(i_x, i_y, std::move(tp));
            }
        }
    }

    // Read the grid
    grid_attach_read_test(get_data(g2, mng_mr));
}

/// This test demonstrates how to call grid buffer without calling host grid
/// object It is especially useful when you don't need to save the objects in
/// host side (e.g. internal spacepoint creation in traccc)
TEST(grids_cuda, grid2_buffer_attach_populator) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    axis2::circular<> xaxis{2u, -1.f, 3.f, mng_mr};
    axis2::regular<> yaxis{2u, 0.f, 6.f, mng_mr};

    grid2_buffer<host_grid2_attach> g2_buffer(
        xaxis, yaxis, {100, 200, 300, 400}, mng_mr, nullptr,
        vecmem::data::buffer_type::resizable);
    copy.setup(g2_buffer._buffer)->wait();

    // Check if the initialization work well
    // Non-zero starting size not working yet so initial argument for sizes is
    // ignored (acts-projects/vecmem#95)
    const auto& ptr = g2_buffer._buffer.host_ptr();
    EXPECT_EQ(ptr[0].size(), 0u);
    EXPECT_EQ(ptr[1].size(), 0u);
    EXPECT_EQ(ptr[2].size(), 0u);
    EXPECT_EQ(ptr[3].size(), 0u);
    EXPECT_EQ(ptr[0].capacity(), 100u);
    EXPECT_EQ(ptr[1].capacity(), 200u);
    EXPECT_EQ(ptr[2].capacity(), 300u);
    EXPECT_EQ(ptr[3].capacity(), 400u);

    // fill each bin with 100 points
    grid_attach_fill_test(g2_buffer);

    host_grid2_attach g2(xaxis, yaxis, mng_mr, point3{0.f, 0.f, 0.f});
    copy(g2_buffer._buffer, g2.data())->wait();

    // Check if each bin has 100 points
    EXPECT_EQ(g2.data()[0].size(), 100u);
    EXPECT_EQ(g2.data()[1].size(), 100u);
    EXPECT_EQ(g2.data()[2].size(), 100u);
    EXPECT_EQ(g2.data()[3].size(), 100u);

    // Check that we can give a non-const buffer to a function expecting
    // a const view.
    grid_attach_read_test(g2_buffer);
}

TEST(grids_cuda, grid2_buffer_attach_populator2) {

    // Helper object for performing memory copies.
    vecmem::cuda::copy copy;

    // memory resource
    vecmem::cuda::managed_memory_resource mng_mr;

    axis2::circular<> xaxis{2u, -1.f, 3.f, mng_mr};
    axis2::regular<> yaxis{2u, 0.f, 6.f, mng_mr};

    grid2_buffer<host_grid2_attach> g2_buffer(xaxis, yaxis, {1, 2, 3, 4},
                                              mng_mr);
    copy.setup(g2_buffer._buffer)->wait();

    // Check if the initialization works well
    const auto& ptr = g2_buffer._buffer.host_ptr();
    EXPECT_EQ(ptr[0].size(), 1u);
    EXPECT_EQ(ptr[1].size(), 2u);
    EXPECT_EQ(ptr[2].size(), 3u);
    EXPECT_EQ(ptr[3].size(), 4u);
    EXPECT_EQ(ptr[0].capacity(), 1u);
    EXPECT_EQ(ptr[1].capacity(), 2u);
    EXPECT_EQ(ptr[2].capacity(), 3u);
    EXPECT_EQ(ptr[3].capacity(), 4u);

    // Assign values to the vector elements
    grid_attach_assign_test(g2_buffer);

    host_grid2_attach g2(xaxis, yaxis, mng_mr, point3{0.f, 0.f, 0.f});
    copy(g2_buffer._buffer, g2.data())->wait();

    // Check the outputs
    auto bin0 = g2.bin(0u);
    EXPECT_EQ(bin0[0], point3({0.f, 1.f, 2.f}));

    auto bin1 = g2.bin(1u);
    EXPECT_EQ(bin1[0], point3({0.f, 1.f, 2.f}));
    EXPECT_EQ(bin1[1], point3({1.f, 2.f, 3.f}));

    auto bin2 = g2.bin(2u);
    EXPECT_EQ(bin2[0], point3({0.f, 1.f, 2.f}));
    EXPECT_EQ(bin2[1], point3({1.f, 2.f, 3.f}));
    EXPECT_EQ(bin2[2], point3({2.f, 3.f, 4.f}));

    auto bin3 = g2.bin(3u);
    EXPECT_EQ(bin3[0], point3({0.f, 1.f, 2.f}));
    EXPECT_EQ(bin3[1], point3({1.f, 2.f, 3.f}));
    EXPECT_EQ(bin3[2], point3({2.f, 3.f, 4.f}));
    EXPECT_EQ(bin3[3], point3({3.f, 4.f, 5.f}));
}
