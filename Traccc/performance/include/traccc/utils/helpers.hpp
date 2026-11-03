// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <string>
#include <string_view>

namespace traccc::plot_helpers {

/// @brief Nested binning struct for booking plots
struct binning {
  /// Constructor with default arguments
  binning(std::string_view b_title = "", int bins = 0, float b_min = 0.f,
          float b_max = 0.f)
      : title(b_title), n_bins(bins), min(b_min), max(b_max) {}

  std::string title;  ///< title to be displayed
  int n_bins;         ///< number of bins
  float min;          ///< minimum value
  float max;          ///< maximum value
};

}  // namespace traccc::plot_helpers
