#pragma once

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

template <typename T,
    unsigned int N_COLS,
    unsigned int N_ROWS >
TRACCC_HOST_DEVICE
bool is_local_maximum(const global_index_t globalIndex,
                      const T value,
                      const vecmem::device_vector<const T>& plane,                      
                      const long long col,
                      const long long row) {
    if (value < 5) return false; // threshold
    static constexpr long long rows = N_ROWS;
    static constexpr long long cols = N_COLS;
     
    for (int dy = -1; dy <= 1; ++dy) {
        const long long neighbor_row = (row + dy + rows) % rows;
        for (int dx = -1; dx <= 1; ++dx) {
            const long long neighbor_col = col + dx;
            if (neighbor_col < 0 || neighbor_col >= cols) continue;

            const auto neighbor_index =
                static_cast<unsigned long long>(neighbor_row) * N_COLS +
                static_cast<unsigned long long>(neighbor_col);

            if (neighbor_index == globalIndex) continue;
            if (plane[neighbor_index] > value) return false;
        } // loop on dx
    } // loop on dy
    return true;
}

} // namespace
