// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <cstdint>
#include <limits>

namespace traccc {

/// Configuration struct for ambiguity resolution
struct ambiguity_resolution_config {
  /// Minimum number of measurement to form a track.
  unsigned int min_meas_per_track = 3;

  /// Max iteration to remove the bad tracks
  unsigned int max_iterations = std::numeric_limits<unsigned int>::max();

  /// Max shared measurements to break the iteration
  unsigned int max_shared_meas = 1;
};

}  // namespace traccc
