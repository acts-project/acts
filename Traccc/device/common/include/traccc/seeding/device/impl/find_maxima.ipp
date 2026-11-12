#pragma once

#include <mutex>

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"
#include "traccc/device/mutex.hpp"
#include "traccc/device/unique_lock.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
void find_maxima(const global_index_t globalIndex,
                 const global_index_t n_wedges,
                 const global_index_t n_cols,
                 const global_index_t n_rows,
                 const vecmem::data::vector_view<const unsigned short>& houghplanes_view,
                 const vecmem::data::vector_view<const unsigned short>& thresholds_view,
                 vecmem::data::jagged_vector_view<global_index_t>& maxima_index_view) {
    const global_index_t dim = n_cols * n_rows;
    if (globalIndex >= dim) return;

    const vecmem::device_vector<const unsigned short> houghplanes(houghplanes_view); 
    const vecmem::device_vector<const unsigned short> thresholds(thresholds_view);
    vecmem::jagged_device_vector<global_index_t> maxima_index(maxima_index_view);

    const auto [col, row] = get_col_row(globalIndex, n_cols, n_rows);
    for (unsigned int wedge_idx(0); wedge_idx<n_wedges; ++wedge_idx) {
        const global_index_t value_index = wedge_idx * dim + globalIndex;
        const unsigned int value = houghplanes[value_index];
        if (value < thresholds[wedge_idx]) continue;

        bool is_maximum = is_local_maximum(globalIndex,
                                           n_cols, n_rows,
                                           value,
                                           &houghplanes[wedge_idx * dim],
                                           col, row);
        if (!is_maximum) continue;
        maxima_index[wedge_idx].push_back(globalIndex);
    } 
}

} // namespace

