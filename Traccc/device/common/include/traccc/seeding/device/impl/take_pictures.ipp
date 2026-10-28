#pragma once

#include <mutex>

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"
#include "traccc/device/mutex.hpp"
#include "traccc/device/unique_lock.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void take_the_picture(const vecmem::jagged_device_vector<unsigned short>::const_reference wedge,
                             vecmem::jagged_device_vector<scalar>::reference image,
                             const global_index_t index,
                             int width, 
                             int height) {
    constexpr int wedge_width = 256;
    constexpr int wedge_height = 7000;

    assert(width > 0 && height > 0);
    
    const auto [col, row] = get_col_row<device::N_COLS, device::N_ROWS>(index);
    const int center_x = static_cast<int>(col);
    const int center_y = static_cast<int>(row);

    for (int ix = 0; ix < width; ++ix) {
        const int x = center_x + ix - width / 2;
        if (x < 0 || x >= wedge_width) continue;

        for (int iy = 0; iy < height; ++iy) {
            int y = (center_y + iy - height / 2) % wedge_height;
            if (y < 0) {
                y += wedge_height;
            }

            const std::size_t output_index =
                static_cast<std::size_t>(ix) * height + iy;
            const std::size_t input_index =
                static_cast<std::size_t>(x) * wedge_height + y;
            image[output_index] = static_cast<scalar>( wedge[input_index] ) / .1;
        }
    }
}


template< typename T,
    unsigned int N_COLS,
    unsigned int N_ROWS >
TRACCC_HOST_DEVICE
void take_pictures(const global_index_t globalIndex,
                   const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::view& houghplanes_view,
                   const vecmem::data::jagged_vector_view<global_index_t>& maxima_index_view,
                   const unsigned int n_maxima,
                   edm::houghplane_collection<scalar, 32, 32>::view& images_view,
                   vecmem::data::vector_view<scalar>& maxima_x_view) {
    if (globalIndex >= n_maxima) return;

    const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device houghplanes(houghplanes_view);
    const vecmem::jagged_device_vector<global_index_t> maxima_index(maxima_index_view);
    edm::houghplane_collection<scalar, 32, 32>::device images(images_view);
    vecmem::device_vector<scalar> maxima_x(maxima_x_view);

    // deduce the correct index
    global_index_t current_size = globalIndex;
    global_index_t wedge_idx = 0;
    global_index_t element_idx = 0;

    for (global_index_t i(0); i<maxima_index.size(); ++i) {
        if (maxima_index[i].size() > current_size) {
           wedge_idx = i;
           element_idx = current_size;
           break;
        }
        current_size -= maxima_index[i].size();
    }
   
    // get variables
    const auto& plane = houghplanes.plane()[wedge_idx];
    const global_index_t maxima_idx = maxima_index[wedge_idx][element_idx];
    take_the_picture(plane, 
                     images.plane()[globalIndex],
                     maxima_idx,
                     32, 32);

    const auto& edges_x = houghplanes.edges_x();
    const auto [col, _] = get_col_row<device::N_COLS, device::N_ROWS>(maxima_idx);
    maxima_x[globalIndex] = 0.50f * (edges_x[col] + edges_x[col + 1]);
}

} // namespace

