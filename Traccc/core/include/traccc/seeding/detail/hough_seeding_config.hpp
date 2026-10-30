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
  scalar min_range_y = -std::numbers::pi_v<scalar>;
  scalar max_range_y = std::numbers::pi_v<scalar>;
  scalar interaction_center = 0.f;
  scalar interaction_window = 150.f;
  scalar eta_width = 0.12f;
  std::array<scalar, 25> wedges = {
  	-2.88f, -2.64f, -2.40f, -2.16f, -1.92f, -1.68f, 
    	-1.44f, -1.20f, -0.96f, -0.72f, -0.48f, -0.24f, 0.00f,
    	0.24f, 0.48f, 0.72f, 0.96f, 1.20f, 1.44f,
    	1.68f, 1.92f, 2.16f, 2.40f, 2.64f, 2.88f
    };
  std::array<unsigned int, 25> min_contributions {
	5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5,
	5, 5, 5, 5, 5, 5
  };
  unsigned long long module_veto = (1ULL << 7) | (1ULL << 8) | (1ULL << 9);
};

}  // namespace traccc
