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
#include "Acts/Utilities/IMultiAxis.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ActsFatras {

std::vector<Segmentizer::ChannelSegment> Segmentizer::segments(
    const Acts::GeometryContext& geoCtx, const Acts::Surface& surface,
    const Acts::IMultiAxis& segmentation, const Segment2D& segment) const {
  // Return if the segmentation is not two-dimensional
  // (strips need to have one bin along the strip)
  if (segmentation.getNAxes() != 2) {
    return {};
  }

  const Acts::IAxis& axis0 = segmentation.getAxis(0);
  const Acts::IAxis& axis1 = segmentation.getAxis(1);

  // Start and end point
  const Acts::Vector2& start = segment[0];
  const Acts::Vector2& end = segment[1];

  // Full path length - the full channel
  auto segment2d = (end - start);
  std::vector<ChannelStep> cSteps;
  Bin2D bstart = {0, 0};
  Bin2D bend = {0, 0};
    
  // Index convention: getBinEdges()[i] == getBinLowerBound(i + 1), and the bin
  // indices here are zero-based (getBin() is one-based, hence the -1 above), so
  // ib <= nBins - 1 and getBinLowerBound(ib + 1) is always in range.
  if (surface.type() == Acts::Surface::SurfaceType::Plane ||
      surface.type() == Acts::Surface::SurfaceType::Cylinder) {
    // For Plane the local frame is Cartesian (x, y); for Cylinder it is the
    // unrolled readout frame (rPhi, z). Either way the cell boundaries are
    // axis-aligned straight lines and the stepping algorithm is identical.
    // Get the segmentation and convert it to lines & arcs
    bstart = {static_cast<unsigned int>(axis0.getBin(start[0]) - 1),
              static_cast<unsigned int>(axis1.getBin(start[1]) - 1)};
    bend = {static_cast<unsigned int>(axis0.getBin(end[0]) - 1),
            static_cast<unsigned int>(axis1.getBin(end[1]) - 1)};
    // Fast single channel exit
    if (bstart == bend) {
      return {ChannelSegment(bstart, {start, end}, segment2d.norm())};
    }
    // The lines channel segment lines along x
    if (bstart[0] != bend[0]) {
      const double k = segment2d.y() / segment2d.x();
      const double d = start.y() - k * start.x();

      const unsigned int xlo = std::min(bstart[0], bend[0]);
      const unsigned int xhi = std::max(bstart[0], bend[0]);
      for (unsigned int ib = xlo + 1; ib <= xhi; ++ib) {
        const double x = axis0.getBinLowerBound(ib + 1);
        cSteps.push_back(ChannelStep{
            {(bstart[0] < bend[0] ? 1 : -1), 0}, {x, k * x + d}, start});
      }
    }
    // The lines channel segment lines along y
    if (bstart[1] != bend[1]) {
      const double k = segment2d.x() / segment2d.y();
      const double d = start.x() - k * start.y();
      const unsigned int ylo = std::min(bstart[1], bend[1]);
      const unsigned int yhi = std::max(bstart[1], bend[1]);
      for (unsigned int ib = ylo + 1; ib <= yhi; ++ib) {
        const double y = axis1.getBinLowerBound(ib + 1);
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
    bstart = {static_cast<unsigned int>(axis0.getBin(pstart[0]) - 1),
              static_cast<unsigned int>(axis1.getBin(pstart[1]) - 1)};
    bend = {static_cast<unsigned int>(axis0.getBin(pend[0]) - 1),
            static_cast<unsigned int>(axis1.getBin(pend[1]) - 1)};

    // Fast single channel exit
    if (bstart == bend) {
      return {ChannelSegment(bstart, {start, end}, segment2d.norm())};
    }

    const double phistart = pstart[1];
    const double phiend = pend[1];

    // The radial boundaries
    if (bstart[0] != bend[0]) {
      const unsigned int rlo = std::min(bstart[0], bend[0]);
      const unsigned int rhi = std::max(bstart[0], bend[0]);
      for (unsigned int ib = rlo + 1; ib <= rhi; ++ib) {
        const double r = axis0.getBinLowerBound(ib + 1);
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
      const unsigned int philo = std::min(bstart[1], bend[1]);
      const unsigned int phihi = std::max(bstart[1], bend[1]);
      for (unsigned int ib = philo + 1; ib <= phihi; ++ib) {
        const double phi = axis1.getBinLowerBound(ib + 1);
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
