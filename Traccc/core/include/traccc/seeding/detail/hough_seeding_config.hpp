// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <numbers>
#include <vector>

#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc {

struct hough_seeding_config {
  unsigned int n_wedges = 25u;
  unsigned int n_cols = 256u;
  unsigned int n_rows = 7000u;
  scalar min_range_x = -1.1f;
  scalar max_range_x = 1.1f;
  scalar min_range_y = -std::numbers::pi_v<scalar>;
  scalar max_range_y = std::numbers::pi_v<scalar>;
  scalar interaction_center = 0.f;
  scalar interaction_window = 150.f;
  scalar eta_width = 0.12f;
  scalar z_window = 550.f;
  std::vector<scalar> wedges = {
      -2.88f, -2.64f, -2.40f, -2.16f, -1.92f, -1.68f, -1.44f, -1.20f, -0.96f,
      -0.72f, -0.48f, -0.24f, 0.00f,  0.24f,  0.48f,  0.72f,  0.96f,  1.20f,
      1.44f,  1.68f,  1.92f,  2.16f,  2.40f,  2.64f,  2.88f};
  std::vector<unsigned short> min_contributions{5, 5, 5, 5, 5, 5, 5, 5, 5,
                                                5, 5, 5, 5, 5, 5, 5, 5, 5,
                                                5, 5, 5, 5, 5, 5, 5};
  unsigned long long module_veto = (1ULL << 7) | (1ULL << 8) | (1ULL << 9);

  static constexpr unsigned int invalid =
      std::numeric_limits<unsigned int>::max();
  std::vector<unsigned int> shift_map = {
      invalid, invalid, invalid,  // 0–2
      0,       invalid, 1,       invalid, 2,       invalid, 3,
      invalid, 4,       invalid, 5,       invalid, 6,            // 3–15
      invalid, invalid,                                          // 16–17
      7,       invalid, 8,       invalid, 9,       invalid, 10,  // 18–24
      invalid, invalid,                                          // 25–26
      11,      invalid, 12,      invalid, 13,      invalid, 14,
      invalid, 15,      invalid, 16,      invalid, 17,  // 27–39
      invalid, invalid, invalid, invalid, invalid,      // 40–44
      18,      invalid, 19,      invalid, 20,      invalid, 21,
      invalid, 22,      invalid, 23,                             // 45–55
      invalid, invalid,                                          // 56–57
      24,      invalid, 25,      invalid, 26,      invalid, 27,  // 58–64
      invalid, invalid,                                          // 65–66
      28,      invalid, 29,      invalid, 30,      invalid, 31,
      invalid, 32,      invalid, 33  // 67–77
  };
};

}  // namespace traccc
