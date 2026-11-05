#pragma once

#include <mutex>

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"
#include "traccc/device/mutex.hpp"
#include "traccc/device/unique_lock.hpp"

namespace traccc::device {

template< typename T>
TRACCC_HOST_DEVICE
void find_maxima(const global_index_t globalIndex,
                 const typename edm::houghplane_collection<T>::view& houghplanes_view,
                 vecmem::data::jagged_vector_view<global_index_t>& maxima_index_view) {
    const global_index_t n_cols = device::N_COLS;
    const global_index_t n_rows = device::N_ROWS;
    if (globalIndex >= n_cols * n_rows) return;

    const typename edm::houghplane_collection<T>::device houghplanes(houghplanes_view);
    vecmem::jagged_device_vector<global_index_t> maxima_index(maxima_index_view);

    const auto [col, row] = get_col_row(globalIndex, n_cols, n_rows);
    for (unsigned int wedge_idx(0); wedge_idx<houghplanes.size(); ++wedge_idx) {
        const auto& plane = houghplanes.plane()[wedge_idx];
        const T value = plane[globalIndex];
        if (value == 0) continue;

        bool is_maximum = is_local_maximum<T>(globalIndex, n_cols, n_rows, value, plane, col, row);
        if (!is_maximum) continue;
        maxima_index[wedge_idx].push_back(globalIndex);
    }    
}

} // namespace

