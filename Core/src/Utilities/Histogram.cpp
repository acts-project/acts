// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Histogram.hpp"

#include <array>
#include <cassert>
#include <string>
#include <vector>

namespace Acts::Experimental {

namespace {

/// Wrap an already projected boost histogram into a Histogram1, carrying over
/// both the axis (including its metadata) and the bin contents.
Histogram1 wrapProjection(const BoostHist& projected, std::string name,
                          std::string title) {
  std::array<AxisVariant, 1> axes = {projected.axis(0)};
  Histogram1 result(std::move(name), std::move(title), axes);

  for (int i = 0; i < projected.axis(0).size(); ++i) {
    result.setBinContent({i}, projected.at(i));
  }

  return result;
}

}  // namespace

// Projection free functions
Histogram1 projectionX(const Histogram2& hist2d) {
  const BoostHist projectedHist = boost::histogram::algorithm::project(
      hist2d.histogram(), std::integral_constant<unsigned, 0>{});

  return wrapProjection(projectedHist, hist2d.name() + "_projX",
                        hist2d.title() + " projection X");
}

Histogram1 projectionY(const Histogram2& hist2d) {
  const BoostHist projectedHist = boost::histogram::algorithm::project(
      hist2d.histogram(), std::integral_constant<unsigned, 1>{});

  return wrapProjection(projectedHist, hist2d.name() + "_projY",
                        hist2d.title() + " projection Y");
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
