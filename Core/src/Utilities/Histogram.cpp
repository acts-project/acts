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

Histogram1 sliceLastAxis(const Histogram2& hist2d, int xBin) {
  const auto& lastAxis = hist2d.histogram().axis(1);
  assert(xBin >= 0 && xBin < hist2d.histogram().axis(0).size() &&
         "x bin index out of range");

  std::array<AxisVariant, 1> axes = {lastAxis};
  Histogram1 slice(hist2d.name() + "_slice_" + std::to_string(xBin),
                   hist2d.title(), axes);

  for (int j = 0; j < lastAxis.size(); ++j) {
    slice.setBinContent({j}, hist2d.binContent({xBin, j}));
  }

  return slice;
}

Histogram1 sliceLastAxis(const Histogram3& hist3d, int xBin, int yBin) {
  const auto& lastAxis = hist3d.histogram().axis(2);
  assert(xBin >= 0 && xBin < hist3d.histogram().axis(0).size() &&
         "x bin index out of range");
  assert(yBin >= 0 && yBin < hist3d.histogram().axis(1).size() &&
         "y bin index out of range");

  std::array<AxisVariant, 1> axes = {lastAxis};
  Histogram1 slice(hist3d.name() + "_slice_" + std::to_string(xBin) + "_" +
                       std::to_string(yBin),
                   hist3d.title(), axes);

  for (int k = 0; k < lastAxis.size(); ++k) {
    slice.setBinContent({k}, hist3d.binContent({xBin, yBin, k}));
  }

  return slice;
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
