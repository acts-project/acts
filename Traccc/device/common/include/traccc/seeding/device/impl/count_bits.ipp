// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/device/global_index.hpp"
#include "traccc/seeding/device/impl/hough_seeding_utilities.ipp"

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void count_bits(
    const global_index_t globalIndex, const global_index_t n_entries,
    const vecmem::data::vector_view<const unsigned long long>& bits_view,
    vecmem::data::vector_view<unsigned short>& entries_view) {
  if (globalIndex >= n_entries) {
    return;
  }

  const vecmem::device_vector<const unsigned long long> bits(bits_view);
  vecmem::device_vector<unsigned short> entries(entries_view);

  const unsigned long long mask = bits[globalIndex];
  if (mask == 0ULL) {
    return;
  }
  entries[globalIndex] = count_bits(mask);
}

}  // namespace traccc::device
