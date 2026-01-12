// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Histogram.hpp"

namespace Acts::Experimental {

// Projection free functions
Histogram1D projectionX(const Histogram2D& hist2d) {
  auto projectedHist = boost::histogram::algorithm::project(
      hist2d.histogram(), std::integral_constant<unsigned, 0>{});

  // Extract single axis from projected histogram
  std::array<AxisVariant, 1> axes = {projectedHist.axis(0)};

  return Histogram1D(hist2d.name() + "_projX", hist2d.title() + " projection X",
                     axes);
}

Histogram1D projectionY(const Histogram2D& hist2d) {
  auto projectedHist = boost::histogram::algorithm::project(
      hist2d.histogram(), std::integral_constant<unsigned, 1>{});

  // Extract single axis from projected histogram
  std::array<AxisVariant, 1> axes = {projectedHist.axis(0)};

  return Histogram1D(hist2d.name() + "_projY", hist2d.title() + " projection Y",
                     axes);
}

}  // namespace Acts::Experimental
