#pragma once

#include <mutex>

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"
#include "traccc/device/mutex.hpp"
#include "traccc/device/unique_lock.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void take_the_picture(const unsigned short* wedge,
                             vecmem::jagged_device_vector<scalar>::reference image,
                             const global_index_t index,
                             const int width, 
                             const int height,
                             const int ix,
                             const int iy) {
    static constexpr scalar normalization = 1.f / 256.f;
    constexpr int wedge_width = 256;
    constexpr int wedge_height = 7000;

    const auto [col, row] = get_col_row(index, device::N_COLS, device::N_ROWS);
    const int center_x = static_cast<int>(col);
    const int center_y = static_cast<int>(row);

    const global_index_t output_index = static_cast<global_index_t>(ix) * height + iy;
    const int x = center_x + ix - width / 2;
    if (x < 0 || x >= wedge_width) {
        image[output_index] = scalar{0};
        return;
    }

    int y = (center_y + iy - height / 2) % wedge_height;
    if (y < 0) y += wedge_height;

    const global_index_t input_index = static_cast<global_index_t>(x) * wedge_height + y;
    image[output_index] =  static_cast<scalar>(wedge[input_index]) * normalization;
}


TRACCC_HOST_DEVICE
void take_pictures(const global_index_t globalIndex,
                   const global_index_t colIndex,
                   const global_index_t rowIndex,
                   const vecmem::data::vector_view<unsigned short>& houghplanes_view,
                   const vecmem::data::vector_view<scalar>& edges_x_view,
                   const vecmem::data::vector_view<scalar>& edges_y_view,
                   const vecmem::data::jagged_vector_view<global_index_t>& maxima_index_view,
                   const unsigned int n_maxima,
                   vecmem::data::jagged_vector_view<scalar>& images_view,
                   vecmem::data::vector_view<scalar>& maxima_x_view) {
    vecmem::device_vector<scalar> edges_x(edges_x_view);
    vecmem::device_vector<scalar> edges_y(edges_y_view);

    const global_index_t n_cols = edges_x.size() - 1;
    const global_index_t n_rows = edges_y.size() - 1;
    
    if (globalIndex >= n_maxima) return;

    const vecmem::device_vector<unsigned short> houghplanes(houghplanes_view);
    const vecmem::jagged_device_vector<global_index_t> maxima_index(maxima_index_view);
    vecmem::jagged_device_vector<scalar> images(images_view);

    // deduce the correct index
    const auto [wedge_idx, element_idx] = get_wedge_element<global_index_t>(globalIndex, maxima_index);    
   
    // get variables
    const global_index_t maxima_idx = maxima_index[wedge_idx][element_idx];
    take_the_picture(&houghplanes[wedge_idx * n_cols * n_rows], 
                     images[globalIndex],
                     maxima_idx,
                     32, 32,
                     colIndex,
                     rowIndex);

    if (colIndex == 0 && rowIndex == 0) {
        vecmem::device_vector<scalar> maxima_x(maxima_x_view);
        const auto [col, _] = get_col_row(maxima_idx, n_cols, n_rows);
        maxima_x[globalIndex] = 0.50f * (edges_x[col] + edges_x[col + 1]);
    }
}

} // namespace

