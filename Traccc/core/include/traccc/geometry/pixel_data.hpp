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
#include "traccc/definitions/qualifiers.hpp"

namespace traccc {

/// A very basic pixel segmentation with
/// a minimum corner and ptich x/y
///
/// No checking on out of bounds done
struct pixel_data {
  scalar min_corner_x = 0.f;
  scalar min_corner_y = 0.f;
  scalar pitch_x = 1.f;
  scalar pitch_y = 1.f;
  char dimension = 2;

  TRACCC_HOST_DEVICE
  vector2 get_pitch() const { return {pitch_x, pitch_y}; };
};

}  // namespace traccc
