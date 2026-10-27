#pragma once

#include <mutex>

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"
#include "traccc/device/mutex.hpp"
#include "traccc/device/unique_lock.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void take_the_picture(const vecmem::device_vector<unsigned short>& wedge,
                             vecmem::device_vector<scalar>& image,
                             const global_index_t col,
                             const global_index_t row,
                             int width, 
                             int height) {
    constexpr int side = 32;
    constexpr int half = side / 2;

    assert(width > 0 && height > 0);
    global_index_t centerX = col;
    global_index_t centerY = row;

    const std::size_t matrixSize =
        static_cast<std::size_t>(width) * height;

    for (int ix = 0; ix < side; ++ix) {
        const int x = centerX + ix - half;

        for (int iy = 0; iy < side; ++iy) {
            int y = (centerY + iy - half) % height;
            if (y < 0) y += height;

            const int outputIndex = ix * side + iy;
            if (x >= 0 && x < width) {
                const std::size_t inputIndex =
                    static_cast<std::size_t>(x) * height + y;

                image[outputIndex] = wedge[inputIndex];
            }
        } // loop on iy
    } // loop on ix
}


template< typename T,
    unsigned int N_COLS,
    unsigned int N_ROWS >
TRACCC_HOST_DEVICE
void take_pictures(const global_index_t globalIndex,
                   const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::view& houghplanes_view,
                   const vecmem::data::jagged_vector_view<global_index_t>& maxima_index_view,
                   const unsigned int n_maxima,
                   edm::houghplane_collection<scalar, 32, 32>::view& images_view) {
    if (globalIndex >= n_maxima) return;

    const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device houghplanes(houghplanes_view);
    const vecmem::jagged_device_vector<global_index_t> maxima_index(maxima_index_view);
    edm::houghplane_collection<scalar, 32, 32>::device images(images_view);
    /*
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
    const auto& plane = houghplanes[wedge_idx];
    const global_index_t maxima_idx = maxima_index[wedge_idx][element_idx];
    const global_index_t maxima_col = maxima_idx / N_COLS;
    const global_index_t maxima_row = maxima_idx % N_COLS;
    */
    // auto target = images[globalIndex];
    // take_the_picture(plane, target, maxima_col, maxima_row, 32, 32);
}

} // namespace

