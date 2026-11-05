#pragma once

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

TRACCC_HOST_DEVICE
void count_maxima(const global_index_t globalIndex,
                  const vecmem::data::vector_view<unsigned short>& houghplanes_view,
                  vecmem::data::vector_view<unsigned int>& number_maxima_view) {
    const unsigned int n_cols = device::N_COLS;
    const unsigned int n_rows = device::N_ROWS;    

    if (globalIndex >= n_cols * n_rows) return;
    const vecmem::device_vector<unsigned short> houghplanes(houghplanes_view);
    vecmem::device_vector<unsigned int> number_maxima(number_maxima_view);

    const auto [col, row] = get_col_row(globalIndex, n_cols, n_rows);
    for (unsigned int wedge_idx(0); wedge_idx<25; ++wedge_idx) {
        const global_index_t value_index = wedge_idx * n_cols * n_rows + globalIndex;
        const unsigned int value = houghplanes[value_index];
        if (value == 0) continue;
        
        bool is_maximum = is_local_maximum(globalIndex, 
                                           n_cols, n_rows, 
                                           value, 
                                           &houghplanes[wedge_idx * n_cols * n_rows],
                                           col, row);
        if (!is_maximum) continue;
        // it is a local maximum
        vecmem::device_atomic_ref<unsigned int> val(number_maxima[wedge_idx]);
        val.fetch_add(1u);
    } // loop on wedges
}

} // namespace
