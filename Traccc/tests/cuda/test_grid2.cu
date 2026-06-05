/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "../../device/cuda/src/utils/cuda_error_handling.hpp"

// Traccc test include(s)
#include "tests/grid2_test.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace traccc {

/*---------------------------------------------------
  test function for grid data with replace populator
  ---------------------------------------------------*/

// cuda kernel for grid_replace_test
__global__ void grid_replace_test_kernel(
    grid2_view<host_grid2_replace> grid_view) {
    // Let's try building the grid object
    device_grid2_replace g2_device(grid_view);

    const auto& axis0 = g2_device.axis_p0();
    const auto& axis1 = g2_device.axis_p1();

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval = (axis1.max - axis1.min) / axis1.n_bins;

    auto gid = threadIdx.x + threadIdx.y * blockDim.x;
    auto pt = point3{axis0.min + gid * x_interval, axis1.min + gid * y_interval,
                     0.5f};

    // replace the bin elements
    g2_device.populate(threadIdx.x, threadIdx.y, std::move(pt));
}

// grid_replace_test implementation
void grid_replace_test(grid2_view<host_grid2_replace> grid_view) {

    const auto& axis0 = grid_view._axis_p0_view;
    const auto& axis1 = grid_view._axis_p1_view;

    unsigned int block_dim = 1;
    dim3 thread_dim(axis0.n_bins, axis1.n_bins);

    // run the kernel
    grid_replace_test_kernel<<<block_dim, thread_dim>>>(grid_view);

    // cuda error check
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// cuda kernel for grid_replace_ci_test
__global__ void grid_replace_ci_test_kernel(
    grid2_view<host_grid2_replace_ci> grid_view) {
    // Let's try building the grid object
    device_grid2_replace_ci g2_device(grid_view);

    const auto& axis0 = g2_device.axis_p0();
    const auto& axis1 = g2_device.axis_p1();

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval =
        axis1.boundaries[threadIdx.y + 1] - axis1.boundaries[threadIdx.y];

    auto gid = threadIdx.x + threadIdx.y * blockDim.x;
    auto pt = point3{axis0.min + gid * x_interval, axis1.min + gid * y_interval,
                     0.5f};

    // replace the bin elements
    g2_device.populate(threadIdx.x, threadIdx.y, std::move(pt));
}

// test function for replace populator with circular and irregular axis
void grid_replace_ci_test(grid2_view<host_grid2_replace_ci> grid_view) {

    const auto& axis0 = grid_view._axis_p0_view;
    const auto& axis1 = grid_view._axis_p1_view;

    unsigned int block_dim = 1;
    dim3 thread_dim(axis0.n_bins, axis1.n_bins);

    // run the kernel
    grid_replace_ci_test_kernel<<<block_dim, thread_dim>>>(grid_view);

    // cuda error check
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/*---------------------------------------------------------------
  test function for grid data with complete populator
  ---------------------------------------------------------------*/

// cuda kernel for grid_complete_test
__global__ void grid_complete_kernel(
    grid2_view<host_grid2_complete> grid_view) {

    // Let's try building the grid object
    device_grid2_complete g2_device(grid_view, point3{0.f, 0.f, 0.f});

    const auto& axis0 = g2_device.axis_p0();
    const auto& axis1 = g2_device.axis_p1();

    auto x_interval = (axis0.max - axis0.min) / axis0.n_bins;
    auto y_interval = (axis1.max - axis1.min) / axis1.n_bins;

    auto bin_id = threadIdx.x + threadIdx.y * blockDim.x;

    for (int i_p = 0; i_p < n_points; i_p++) {
        auto gid = i_p + bin_id * n_points;
        auto pt = point3{axis0.min + gid * x_interval,
                         axis1.min + gid * y_interval, 0.5f};
        g2_device.populate(threadIdx.x, threadIdx.y, std::move(pt));
    }
}

// grid_complete_test implementation
void grid_complete_test(grid2_view<host_grid2_complete> grid_view) {

    const auto& axis0 = grid_view._axis_p0_view;
    const auto& axis1 = grid_view._axis_p1_view;

    unsigned int block_dim = 1;
    dim3 thread_dim(axis0.n_bins, axis1.n_bins);

    // run the kernel
    grid_complete_kernel<<<block_dim, thread_dim>>>(grid_view);

    // cuda error check
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/*---------------------------------------------------------
  read test function for grid with attach populator
  ---------------------------------------------------------*/

// cuda kernel for attach_read_test
__global__ void grid_attach_read_test_kernel(
    const_grid2_view<host_grid2_attach> grid_view) {

    // Let's try building the grid object
    const const_device_grid2_attach g2_device(grid_view, point3{0.f, 0.f, 0.f});

    auto data = g2_device.bin(threadIdx.x, threadIdx.y);

    for (auto& pt : data) {
        // printf("%f %f %f \n", pt[0], pt[1], pt[2]);
    }
}

// attach_read_test implementation
void grid_attach_read_test(const_grid2_view<host_grid2_attach> grid_view) {

    const auto& axis0 = grid_view._axis_p0_view;
    const auto& axis1 = grid_view._axis_p1_view;

    unsigned int block_dim = 1;
    dim3 thread_dim(axis0.n_bins, axis1.n_bins);

    // run the kernel
    grid_attach_read_test_kernel<<<block_dim, thread_dim>>>(grid_view);

    // cuda error check
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/*---------------------------------------------------------
  fill test function for grid buffer with attach populator
  ---------------------------------------------------------*/

// cuda kernel for attach_fill_test
__global__ void grid_attach_fill_test_kernel(
    grid2_view<host_grid2_attach> grid_view) {

    // Let's try building the grid object
    device_grid2_attach g2_device(grid_view);

    // Fill with 100 points
    auto pt = point3{scalar(1.) * threadIdx.x, scalar(1.) * threadIdx.x,
                     scalar(1.) * threadIdx.x};
    g2_device.populate(blockIdx.x, blockIdx.y, std::move(pt));

    __syncthreads();

    if (threadIdx.x == 0 && blockIdx.x == 0 && blockIdx.y == 0) {
        auto pts = g2_device.bin(blockIdx.x, blockIdx.y);
        for (int i = 0; i < 100; i++) {
            // printf("%f %f %f \n", pts[i][0], pts[i][1], pts[i][2]);
        }
    }
}

// attach_fill_test implementation
void grid_attach_fill_test(grid2_view<host_grid2_attach> grid_view) {

    const auto& axis0 = grid_view._axis_p0_view;
    const auto& axis1 = grid_view._axis_p1_view;

    dim3 block_dim(axis0.n_bins, axis1.n_bins);
    unsigned int thread_dim = 100;

    // run the kernel
    grid_attach_fill_test_kernel<<<block_dim, thread_dim>>>(grid_view);

    // cuda error check
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

/*---------------------------------------------------------
  assign test function for grid buffer with attach populator
  ---------------------------------------------------------*/

__global__ void grid_attach_assign_test_kernel(
    grid2_view<host_grid2_attach> grid_view) {

    // Let's try building the grid object
    device_grid2_attach g2_device(grid_view);

    auto pts = g2_device.bin(threadIdx.x, threadIdx.y);

    for (std::size_t i = 0u; i < pts.size(); i++) {
        pts[i] = {static_cast<scalar>(i), static_cast<scalar>(i + 1u),
                  static_cast<scalar>(i + 2u)};
    }
}

// attach_fill_test implementation
void grid_attach_assign_test(grid2_view<host_grid2_attach> grid_view) {

    unsigned int block_dim = 1;
    dim3 thread_dim(2, 2);

    // run the kernel
    grid_attach_assign_test_kernel<<<block_dim, thread_dim>>>(grid_view);

    // cuda error check
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
    TRACCC_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace traccc
