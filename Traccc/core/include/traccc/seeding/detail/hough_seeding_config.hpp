// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <numbers>

#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc {

struct hough_seeding_config {
  scalar min_range_x = -1.1f;
  scalar max_range_x = 1.1f;
  scalar min_range_y = -std::numbers::pi_v<float>;
  scalar max_range_y = std::numbers::pi_v<float>;
  scalar interaction_center = 0.f;
  scalar interaction_window = 150.f;
};

}  // namespace traccc
