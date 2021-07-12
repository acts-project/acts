// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"

#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

Acts::TGeoCylinderDiscSplitter::TGeoCylinderDiscSplitter(
    const TGeoCylinderDiscSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoCylinderDiscSplitter::split(
    const GeometryContext& gctx,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  // Segments are detected, attempt a split
  if (m_cfg.circularSegments > 0 or m_cfg.regularSegments > 0) {
    const Acts::Surface& sf = tgde->surface();

    // Thickness
    auto tgIdentifier = tgde->identifier();
    const TGeoNode& tgNode = tgde->tgeoNode();
    ActsScalar tgThickness = tgde->thickness();

    // Splitting for discs
    if (sf.type() == Acts::Surface::Disc and
        sf.bounds().type() == Acts::SurfaceBounds::eDisc) {
      std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
          tgDetectorElements = {};
      tgDetectorElements.reserve(std::abs(m_cfg.circularSegments) *
                                 std::abs(m_cfg.regularSegments));

      const Acts::Vector3 discCenter = sf.center(gctx);

      const auto& boundValues = sf.bounds().values();
      ActsScalar discMinR = boundValues[Acts::RadialBounds::eMinR];
      ActsScalar discMaxR = boundValues[Acts::RadialBounds::eMaxR];

      ActsScalar phiStep = 2 * M_PI / m_cfg.circularSegments;
      ActsScalar cosPhiHalf = std::cos(0.5 * phiStep);
      ActsScalar sinPhiHalf = std::sin(0.5 * phiStep);

      std::vector<ActsScalar> radialValues = {};
      if (m_cfg.regularSegments > 1) {
        ActsScalar rStep = (discMaxR - discMinR) / m_cfg.regularSegments;
        radialValues.reserve(m_cfg.regularSegments);
        for (int ir = 0; ir <= m_cfg.regularSegments; ++ir) {
          radialValues.push_back(discMinR + ir * rStep);
        }
      } else {
        radialValues = {discMinR, discMaxR};
      }

      for (size_t ir = 1; ir < radialValues.size(); ++ir) {
        ActsScalar minR = radialValues[ir - 1];
        ActsScalar maxR = radialValues[ir];

        ActsScalar maxLocY = maxR * cosPhiHalf;
        ActsScalar minLocY = minR * cosPhiHalf;

        ActsScalar hR = 0.5 * (maxLocY + minLocY);
        ActsScalar hY = 0.5 * (maxLocY - minLocY);
        ActsScalar hXminY = minR * sinPhiHalf;
        ActsScalar hXmaxY = maxR * sinPhiHalf;

        auto tgTrapezoid =
            std::make_shared<Acts::TrapezoidBounds>(hXminY, hXmaxY, hY);

        for (int im = 0; im < m_cfg.circularSegments; ++im) {
          // Get the moduleTransform
          double phi = -M_PI + im * phiStep;
          auto tgTransform =
              Transform3(Translation3(hR * std::cos(phi), hR * std::sin(phi),
                                      discCenter.z()) *
                         AngleAxis3(phi - 0.5 * M_PI, Vector3::UnitZ()));

          // Create a new detector element per split
          auto tgDetectorElement = std::make_shared<Acts::TGeoDetectorElement>(
              tgIdentifier, tgNode, tgTransform, tgTrapezoid, tgThickness);

          tgDetectorElements.push_back(tgDetectorElement);
        }
      }

      return tgDetectorElements;
    }
  }
  return {tgde};
}
