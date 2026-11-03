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

/// Payload for the @c device::update_tip_length_buffer function
struct update_tip_length_buffer_payload {
  vecmem::data::vector_view<const unsigned int> old_tip_length;
  vecmem::data::vector_view<unsigned int> new_tip_length;
  vecmem::data::vector_view<const unsigned int> measurement_votes;
  vecmem::data::vector_view<unsigned int> tip_to_output_map;
  float min_measurement_voting_fraction;
};

TRACCC_HOST_DEVICE inline void update_tip_length_buffer(
    global_index_t thread_id, const update_tip_length_buffer_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/finding/device/impl/update_tip_length_buffer.ipp"
