#pragma once

#include "traccc/device/global_index.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE
void count_bits(const global_index_t globalIndex,
                const unsigned long long n_entries,
                const vecmem::data::vector_view<unsigned long long>& bits_view,
                vecmem::data::vector_view<unsigned short>& entries_view) {
    if (globalIndex >= n_entries) return;

    const vecmem::device_vector<unsigned long long> bits(bits_view);
    vecmem::device_vector<unsigned short> entries(entries_view);
   
    unsigned long long mask = bits[globalIndex]; 
    if (mask == 0ULL) return;
    entries[globalIndex] = count_bits(mask);
}

} // namespace

