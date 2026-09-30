// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <cstddef>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::sort_tracks_per_measurement
/// function
struct sort_tracks_per_measurement_payload {
  /**
   * @brief View object to the tracks per measurement
   */
  vecmem::data::jagged_vector_view<unsigned int> tracks_per_measurement_view;
};

}  // namespace traccc::device
