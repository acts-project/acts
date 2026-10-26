#pragma once

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

template< typename T,
    unsigned int N_COLS,
    unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void find_maxima(const global_index_t globalIndex,
                 const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::view& houghplanes_view,
                 vecmem::data::jagged_vector_view<unsigned long long>& maxima_index_view,
                 vecmem::data::jagged_vector_view<scalar>& maxima_x_values_view,
                 edm::houghplane_collection<scalar, 32, 32>::view& images_view) {
    if (globalIndex >= N_COLS * N_ROWS) return;
    /*
    const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device houghplanes(houghplanes_view);
    vecmem::device_vector<unsigned int> number_maxima(number_maxima_view);

    const long long col = globalIndex % N_COLS;
    const long long row = globalIndex / N_COLS;
    constexpr long long rows = N_ROWS;
    constexpr long long cols = N_COLS;

    for (unsigned int wedge_idx(0); wedge_idx<houghplanes.size(); ++wedge_idx) {
        const auto& plane = houghplanes.plane()[wedge_idx];
        const T value = plane[globalIndex];
        if (value == 0) continue;
        
        bool is_maximum = true;
        for (int dy = -1; dy <= 1; ++dy) {
            if (!is_maximum) break;
            const long long neighbor_row = (row + dy + rows) % rows;

            for (int dx = -1; dx <= 1; ++dx) {
                const long long neighbor_col = col + dx;
                if (neighbor_col < 0 || neighbor_col >= cols) continue;

                const auto neighbor_index =
                    static_cast<unsigned long long>(neighbor_row) * N_COLS +
                    static_cast<unsigned long long>(neighbor_col);

                if (neighbor_index == globalIndex) continue;
                if (plane[neighbor_index] > value) {
                    is_maximum = false;
                    break;
                }
            } // loop on dx
        } // loop on dy

        // it is a local maximum
        if (!is_maximum) continue;
        vecmem::device_atomic_ref<unsigned int> val(number_maxima[wedge_idx]);
        val.fetch_add(1u);
    } // loop on wedges
    */
}

} // namespace
