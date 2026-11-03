// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Payload structure for the @c device::gather_measurement_votes function
struct gather_measurement_votes_payload {
  vecmem::data::vector_view<const unsigned long long int> insertion_mutex;
  vecmem::data::vector_view<const unsigned int> tip_index;
  vecmem::data::vector_view<unsigned int> votes_per_tip;
  unsigned int max_num_tracks_per_measurement;
};

TRACCC_HOST_DEVICE inline void gather_measurement_votes(
    global_index_t thread_id, const gather_measurement_votes_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/gather_measurement_votes.ipp"
