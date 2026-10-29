#pragma once

#include "traccc/device/global_index.hpp"
#include "vecmem/memory/device_atomic_ref.hpp"


namespace traccc::device {

template <unsigned int N_COLS,
          unsigned int N_ROWS>
TRACCC_HOST_DEVICE
inline std::pair<unsigned int, unsigned int> get_col_row(unsigned long long global_index) {
    assert(global_index < N_COLS * N_ROWS);
    return {static_cast<unsigned int>(global_index / N_ROWS),
            static_cast<unsigned int>(global_index % N_ROWS)};
}

TRACCC_HOST_DEVICE
inline unsigned short count_bits(unsigned long long x) {
    unsigned short nOnes = 0;
    for (int i(0); i<64; ++i) {
       if ((x & (1ULL << i)) != 0ULL) ++nOnes;
    }
    return nOnes;
}

TRACCC_HOST_DEVICE
int quant(scalar min, scalar max, unsigned nSteps, scalar val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}

TRACCC_HOST_DEVICE
scalar transform(const scalar qopt_value,
                 const scalar r,
                 const scalar phi,
         const scalar magnetic_field_z) {
    const scalar d0 = 0;
    scalar phi_value = asin(0.50f * r * magnetic_field_z * qopt_value - d0 / r) + phi;
    ::remainderf(phi_value, 2.00f * std::numbers::pi_v<scalar>);
    return phi_value;
}

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

    for (int dx = -1; dx <= 1; ++dx) {
            const long long neighbor_col = col + dx;
            if (neighbor_col < 0 || neighbor_col >= cols) continue;     
        for (int dy = -1; dy <= 1; ++dy) {
            const long long neighbor_row = (row + dy + rows) % rows;

            const auto neighbor_index =
                static_cast<unsigned long long>(neighbor_col) * N_ROWS +
                static_cast<unsigned long long>(neighbor_row);

            if (neighbor_index == globalIndex) continue;
            if (plane[neighbor_index] > value) return false;
        } // loop on dx
    } // loop on dy
    return true;
}

} // namespace
