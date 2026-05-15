// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Histogram.hpp"

#include <cassert>
#include <vector>

namespace Acts::Experimental {

// Projection free functions
Histogram1 projectionX(const Histogram2& hist2d) {
  auto projectedHist = boost::histogram::algorithm::project(
      hist2d.histogram(), std::integral_constant<unsigned, 0>{});

  // Extract single axis from projected histogram
  std::array<AxisVariant, 1> axes = {projectedHist.axis(0)};

  return Histogram1(hist2d.name() + "_projX", hist2d.title() + " projection X",
                    axes);
}

Histogram1 projectionY(const Histogram2& hist2d) {
  auto projectedHist = boost::histogram::algorithm::project(
      hist2d.histogram(), std::integral_constant<unsigned, 1>{});

  // Extract single axis from projected histogram
  std::array<AxisVariant, 1> axes = {projectedHist.axis(0)};

  return Histogram1(hist2d.name() + "_projY", hist2d.title() + " projection Y",
                    axes);
}

std::vector<double> extractBinEdges(const AxisVariant& axis) {
  assert(axis.size() > 0 && "Axis must have at least one bin");
  std::vector<double> edges(axis.size() + 1);
  for (int i = 0; i < axis.size(); ++i) {
    edges.at(i) = axis.bin(i).lower();
  }
  edges.back() = axis.bin(axis.size() - 1).upper();

  return edges;
}

}  // namespace Acts::Experimental
