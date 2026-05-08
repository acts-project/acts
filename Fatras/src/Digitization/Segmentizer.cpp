// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/Segmentizer.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/IntersectionHelper2D.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

namespace ActsFatras {

std::vector<Segmentizer::ChannelSegment> Segmentizer::segments(
    const Acts::GeometryContext& geoCtx, const Acts::Surface& surface,
    const std::vector<Acts::DirectedProtoAxis>& segmentation,
    const Segment2D& segment) const {
  // Return if the segmentation is not two-dimensional
  // (strips need to have one bin along the strip)
  if (segmentation.size() != 2) {
    return {};
  }

  // Start and end point
  const Acts::Vector2& start = segment[0];
  const Acts::Vector2& end = segment[1];

  // Full path length - the full channel
  auto segment2d = (end - start);
  std::vector<ChannelStep> cSteps;
  Bin2D bstart = {0, 0};
  Bin2D bend = {0, 0};

  if (surface.type() == Acts::Surface::SurfaceType::Plane ||
      surface.type() == Acts::Surface::SurfaceType::Cylinder) {
    // For Plane the local frame is Cartesian (x, y); for Cylinder it is the
    // unrolled readout frame (rPhi, z). Either way the cell boundaries are
    // axis-aligned straight lines and the stepping algorithm is identical.
    // Get the segmentation and convert it to lines & arcs
    bstart = {static_cast<unsigned int>(segmentation.at(0).bin(start[0]) - 1),
              static_cast<unsigned int>(segmentation.at(1).bin(start[1]) - 1)};
    bend = {static_cast<unsigned int>(segmentation.at(0).bin(end[0]) - 1),
            static_cast<unsigned int>(segmentation.at(1).bin(end[1]) - 1)};
    // Fast single channel exit
    if (bstart == bend) {
      return {ChannelSegment(bstart, {start, end}, segment2d.norm())};
    }
    // The lines channel segment lines along x
    if (bstart[0] != bend[0]) {
      const double k = segment2d.y() / segment2d.x();
      const double d = start.y() - k * start.x();

      const std::vector<double> xboundaries = segmentation.at(0).binEdges();
      const std::span<const double> xbbounds(
          xboundaries.begin() + std::min(bstart[0], bend[0]) + 1,
          xboundaries.begin() + std::max(bstart[0], bend[0]) + 1);
      for (const double x : xbbounds) {
        cSteps.push_back(ChannelStep{
            {(bstart[0] < bend[0] ? 1 : -1), 0}, {x, k * x + d}, start});
      }
    }
    // The lines channel segment lines along y
    if (bstart[1] != bend[1]) {
      const double k = segment2d.x() / segment2d.y();
      const double d = start.x() - k * start.y();
      const std::vector<double> yboundaries = segmentation.at(1).binEdges();
      const std::span<const double> ybbounds(
          yboundaries.begin() + std::min(bstart[1], bend[1]) + 1,
          yboundaries.begin() + std::max(bstart[1], bend[1]) + 1);
      for (const double y : ybbounds) {
        cSteps.push_back(ChannelStep{
            {0, (bstart[1] < bend[1] ? 1 : -1)}, {k * y + d, y}, start});
      }
    }

  } else if (surface.type() == Acts::Surface::SurfaceType::Disc) {
    const Acts::Vector2 pstart(Acts::VectorHelpers::perp(start),
                               Acts::VectorHelpers::phi(start));
    const Acts::Vector2 pend(Acts::VectorHelpers::perp(end),
                             Acts::VectorHelpers::phi(end));

    // Get the segmentation and convert it to lines & arcs
    bstart = {static_cast<unsigned int>(segmentation.at(0).bin(pstart[0]) - 1),
              static_cast<unsigned int>(segmentation.at(1).bin(pstart[1]) - 1)};
    bend = {static_cast<unsigned int>(segmentation.at(0).bin(pend[0]) - 1),
            static_cast<unsigned int>(segmentation.at(1).bin(pend[1]) - 1)};

    // Fast single channel exit
    if (bstart == bend) {
      return {ChannelSegment(bstart, {start, end}, segment2d.norm())};
    }

    const double phistart = pstart[1];
    const double phiend = pend[1];

    // The radial boundaries
    if (bstart[0] != bend[0]) {
      const std::vector<double> rboundaries = segmentation.at(0).binEdges();
      const std::span<const double> rbbounds(
          rboundaries.begin() + std::min(bstart[0], bend[0]) + 1,
          rboundaries.begin() + std::max(bstart[0], bend[0]) + 1);
      for (const double r : rbbounds) {
        const auto radIntersection =
            Acts::detail::IntersectionHelper2D::intersectCircleSegment(
                r, std::min(phistart, phiend), std::max(phistart, phiend),
                start, (end - start).normalized());
        cSteps.push_back(ChannelStep{{(bstart[0] < bend[0] ? 1 : -1), 0},
                                     radIntersection.position(),
                                     start});
      }
    }
    // The phi boundaries
    if (bstart[1] != bend[1]) {
      const double referenceR =
          surface.referencePositionValue(geoCtx, Acts::AxisDirection::AxisR);
      const Acts::Vector2 origin = {0., 0.};
      const std::vector<double> phiboundaries = segmentation.at(1).binEdges();
      const std::span<const double> phibbounds(
          phiboundaries.begin() + std::min(bstart[1], bend[1]) + 1,
          phiboundaries.begin() + std::max(bstart[1], bend[1]) + 1);

      for (const double phi : phibbounds) {
        Acts::Vector2 philine(referenceR * std::cos(phi),
                              referenceR * std::sin(phi));
        const auto phiIntersection =
            Acts::detail::IntersectionHelper2D::intersectSegment(
                origin, philine, start, (end - start).normalized());
        cSteps.push_back(ChannelStep{{0, (bstart[1] < bend[1] ? 1 : -1)},
                                     phiIntersection.position(),
                                     start});
      }
    }
  }

  // Register the last step if successful
  if (!cSteps.empty()) {
    cSteps.push_back(ChannelStep({0, 0}, end, start));
    std::ranges::sort(cSteps, std::less<ChannelStep>{});
  }

  std::vector<ChannelSegment> cSegments;
  cSegments.reserve(cSteps.size());

  Bin2D currentBin = {bstart[0], bstart[1]};
  BinDelta2D lastDelta = {0, 0};
  Acts::Vector2 lastIntersect = start;
  double lastPath = 0.;
  for (auto& cStep : cSteps) {
    currentBin[0] += lastDelta[0];
    currentBin[1] += lastDelta[1];
    double path = cStep.path - lastPath;
    cSegments.push_back(
        ChannelSegment(currentBin, {lastIntersect, cStep.intersect}, path));
    lastPath = cStep.path;
    lastDelta = cStep.delta;
    lastIntersect = cStep.intersect;
  }

  return cSegments;
}

}  // namespace ActsFatras
