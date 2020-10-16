// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/Channelizer.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"
#include <Acts/Utilities/BinUtility.hpp>

std::vector<ActsFatras::Channelizer::ChannelSegment>
ActsFatras::Channelizer::segments(const ActsFatras::DigitizationInput& dInput,
                                  const Acts::Vector2D& start,
                                  const Acts::Vector2D& end) const {
  // Return if the segmentation is not two-dimensional
  // (strips need to have one bin along the strip)
  const auto& segmentation = dInput.segmentation;
  if (segmentation.dimensions() != 2) {
    return {};
  }

  // Full path length - the full channel
  auto segment2d = (end - start);

  std::vector<ChannelStep> cSteps;
  unsigned int b0 = 0;
  unsigned int e0 = 0;
  unsigned int b1 = 0;
  unsigned int e1 = 0;

  // The surface properties
  const auto& surface = dInput.surface;
  if (surface->type() == Acts::Surface::SurfaceType::Plane) {
    // Get the segmentation and convert it to lines & arcs
    b0 = segmentation.bin(start, 0);
    e0 = segmentation.bin(end, 0);
    b1 = segmentation.bin(start, 1);
    e1 = segmentation.bin(end, 1);
    // Fast single channel exit
    if (b0 == e0 and b1 == e1) {
      return {ChannelSegment({b0, b1}, segment2d.norm())};
    }
    // The lines channel segment lines along x
    if (b0 != e0) {
      double k = segment2d.y() / segment2d.x();
      double d = start.y() - k * start.x();

      const auto& xboundaries = segmentation.binningData()[0].boundaries();
      std::vector<double> xbbounds = {
          xboundaries.begin() + std::min(b0, e0) + 1,
          xboundaries.begin() + std::max(b0, e0) + 1};
      for (const auto x : xbbounds) {
        cSteps.push_back(
            ChannelStep{{(b0 < e0 ? 1 : -1), 0}, {x, k * x + d}, start});
      }
    }
    // The lines channel segment lines along y
    if (b1 != e1) {
      double k = segment2d.x() / segment2d.y();
      double d = start.x() - k * start.y();
      const auto& yboundaries = segmentation.binningData()[1].boundaries();
      std::vector<double> ybbounds = {
          yboundaries.begin() + std::min(b1, e1) + 1,
          yboundaries.begin() + std::max(b1, e1) + 1};
      for (const auto y : ybbounds) {
        cSteps.push_back(
            ChannelStep{{0, (b1 < e1 ? 1 : -1)}, {k * y + d, y}, start});
      }
    }

  } else if (surface->type() == Acts::Surface::SurfaceType::Disc) {
  }

  // Register the last step
  cSteps.push_back(ChannelStep({0, 0}, end, start));
  std::sort(cSteps.begin(), cSteps.end());

  std::vector<ChannelSegment> cSegments;
  cSegments.reserve(cSteps.size());

  std::array<unsigned int, 2> cbin = {b0, b1};
  std::array<int, 2> cdelta = {0, 0};
  double cpath = 0.;
  for (auto& cStep : cSteps) {
    cbin[0] += cdelta[0];
    cbin[1] += cdelta[1];
    double path = cStep.path - cpath;
    cSegments.push_back(ChannelSegment(cbin, path));
    cpath = cStep.path;
    cdelta = cStep.delta;
  }

  return cSegments;
}
