#pragma once

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

template< typename T,
    unsigned int N_COLS,
    unsigned int N_ROWS>
TRACCC_HOST_DEVICE
void count_maxima(const global_index_t globalIndex,
                  const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::view& houghplanes_view,
                  vecmem::data::vector_view<unsigned int>& number_maxima_view) {
    if (globalIndex >= N_COLS * N_ROWS) return;
    const typename edm::houghplane_collection<T, N_COLS, N_ROWS>::device houghplanes(houghplanes_view);
    vecmem::device_vector<unsigned int> number_maxima(number_maxima_view);

    const long long col = globalIndex % N_COLS;
    const long long row = globalIndex / N_COLS;

    for (unsigned int wedge_idx(0); wedge_idx<houghplanes.size(); ++wedge_idx) {
        const auto& plane = houghplanes.plane()[wedge_idx];
        const T value = plane[globalIndex];
        if (value == 0) continue;
        
        bool is_maximum = is_local_maximum<T, N_COLS, N_ROWS>(globalIndex, value, plane, col, row);
        if (!is_maximum) continue;
        // it is a local maximum
        vecmem::device_atomic_ref<unsigned int> val(number_maxima[wedge_idx]);
        val.fetch_add(1u);
    } // loop on wedges

}

} // namespace
