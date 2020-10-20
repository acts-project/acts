// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/Channelizer.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/IntersectionHelper2D.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"

std::vector<ActsFatras::Channelizer::ChannelSegment>
ActsFatras::Channelizer::segments(const Acts::GeometryContext& geoCtx,
                                  const Acts::Surface& surface,
                                  const Acts::BinUtility& segmentation,
                                  const Acts::Vector2D& start,
                                  const Acts::Vector2D& end) const {
  // Return if the segmentation is not two-dimensional
  // (strips need to have one bin along the strip)
  if (segmentation.dimensions() != 2) {
    return {};
  }

  // Full path length - the full channel
  auto segment2d = (end - start);
  std::vector<ChannelStep> cSteps;
  std::array<unsigned int, 2> bstart = {0, 0};
  std::array<unsigned int, 2> bend = {0, 0};

  if (surface.type() == Acts::Surface::SurfaceType::Plane) {
    // Get the segmentation and convert it to lines & arcs
    bstart = {static_cast<unsigned int>(segmentation.bin(start, 0)),
              static_cast<unsigned int>(segmentation.bin(start, 1))};
    bend = {static_cast<unsigned int>(segmentation.bin(end, 0)),
            static_cast<unsigned int>(segmentation.bin(end, 1))};
    // Fast single channel exit
    if (bstart == bend) {
      return {ChannelSegment(bstart, segment2d.norm())};
    }
    // The lines channel segment lines along x
    if (bstart[0] != bend[0]) {
      double k = segment2d.y() / segment2d.x();
      double d = start.y() - k * start.x();

      const auto& xboundaries = segmentation.binningData()[0].boundaries();
      std::vector<double> xbbounds = {
          xboundaries.begin() + std::min(bstart[0], bend[0]) + 1,
          xboundaries.begin() + std::max(bstart[0], bend[0]) + 1};
      for (const auto x : xbbounds) {
        cSteps.push_back(ChannelStep{
            {(bstart[0] < bend[0] ? 1 : -1), 0}, {x, k * x + d}, start});
      }
    }
    // The lines channel segment lines along y
    if (bstart[1] != bend[1]) {
      double k = segment2d.x() / segment2d.y();
      double d = start.x() - k * start.y();
      const auto& yboundaries = segmentation.binningData()[1].boundaries();
      std::vector<double> ybbounds = {
          yboundaries.begin() + std::min(bstart[1], bend[1]) + 1,
          yboundaries.begin() + std::max(bstart[1], bend[1]) + 1};
      for (const auto y : ybbounds) {
        cSteps.push_back(ChannelStep{
            {0, (bstart[1] < bend[1] ? 1 : -1)}, {k * y + d, y}, start});
      }
    }

  } else if (surface.type() == Acts::Surface::SurfaceType::Disc) {
    Acts::Vector2D pstart(Acts::VectorHelpers::perp(start),
                          Acts::VectorHelpers::phi(start));
    Acts::Vector2D pend(Acts::VectorHelpers::perp(end),
                        Acts::VectorHelpers::phi(end));

    // Get the segmentation and convert it to lines & arcs
    bstart = {static_cast<unsigned int>(segmentation.bin(pstart, 0)),
              static_cast<unsigned int>(segmentation.bin(pstart, 1))};
    bend = {static_cast<unsigned int>(segmentation.bin(pend, 0)),
            static_cast<unsigned int>(segmentation.bin(pend, 1))};

    // Fast single channel exit
    if (bstart == bend) {
      return {ChannelSegment(bstart, segment2d.norm())};
    }

    double phistart = pstart[1];
    double phiend = pend[1];

    // The radial boundaries
    if (bstart[0] != bend[0]) {
      const auto& rboundaries = segmentation.binningData()[0].boundaries();
      std::vector<double> rbbounds = {
          rboundaries.begin() + std::min(bstart[0], bend[0]) + 1,
          rboundaries.begin() + std::max(bstart[0], bend[0]) + 1};
      for (const auto& r : rbbounds) {
        auto radIntersection =
            Acts::detail::IntersectionHelper2D::intersectCircleSegment(
                r, std::min(phistart, phiend), std::max(phistart, phiend),
                start, (end - start).normalized());
        cSteps.push_back(ChannelStep{{(bstart[0] < bend[0] ? 1 : -1), 0},
                                     radIntersection.position,
                                     start});
      }
    }
    // The phi boundaries
    if (bstart[1] != bend[1]) {
      double referenceR = surface.binningPositionValue(geoCtx, Acts::binR);
      Acts::Vector2D origin = {0., 0.};
      const auto& phiboundaries = segmentation.binningData()[1].boundaries();
      std::vector<double> phibbounds = {
          phiboundaries.begin() + std::min(bstart[1], bend[1]) + 1,
          phiboundaries.begin() + std::max(bstart[1], bend[1]) + 1};

      for (const auto& phi : phibbounds) {
        Acts::Vector2D philine(referenceR * std::cos(phi),
                               referenceR * std::sin(phi));
        auto phiIntersection =
            Acts::detail::IntersectionHelper2D::intersectSegment(
                origin, philine, start, (end - start).normalized());
        cSteps.push_back(ChannelStep{{0, (bstart[1] < bend[1] ? 1 : -1)},
                                     phiIntersection.position,
                                     start});
      }
    }
  }

  // Register the last step if successful
  if (not cSteps.empty()) {
    cSteps.push_back(ChannelStep({0, 0}, end, start));
    std::sort(cSteps.begin(), cSteps.end());
  }

  std::vector<ChannelSegment> cSegments;
  cSegments.reserve(cSteps.size());

  std::array<unsigned int, 2> cbin = {bstart[0], bstart[1]};
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
