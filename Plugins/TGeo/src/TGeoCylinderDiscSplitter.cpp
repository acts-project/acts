// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cmath>
#include <cstdlib>
#include <numbers>
#include <utility>

Acts::TGeoCylinderDiscSplitter::TGeoCylinderDiscSplitter(
    const TGeoCylinderDiscSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoCylinderDiscSplitter::split(
    const GeometryContext& gctx,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  const Acts::Surface& sf = tgde->surface();
  // Thickness
  auto tgIdentifier = tgde->identifier();
  const TGeoNode& tgNode = tgde->tgeoNode();
  long double tgThickness = tgde->thickness();

  // Disc segments are detected, attempt a split
  if (m_cfg.discPhiSegments > 0 || m_cfg.discRadialSegments > 0) {
    // Splitting for discs detected
    if (sf.type() == Acts::Surface::Disc &&
        sf.bounds().type() == Acts::SurfaceBounds::eDisc) {
      ACTS_DEBUG("- splitting detected for a Disc shaped sensor.");

      std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
          tgDetectorElements = {};
      tgDetectorElements.reserve(std::abs(m_cfg.discPhiSegments) *
                                 std::abs(m_cfg.discRadialSegments));

      const Acts::Vector3 discCenter = sf.center(gctx);

      const auto& boundValues = sf.bounds().values();
      long double discMinR = boundValues[Acts::RadialBounds::eMinR];
      long double discMaxR = boundValues[Acts::RadialBounds::eMaxR];

      long double phiStep = 2 * std::numbers::pi / m_cfg.discPhiSegments;
      long double cosPhiHalf = std::cos(0.5 * phiStep);
      long double sinPhiHalf = std::sin(0.5 * phiStep);

      std::vector<long double> radialValues = {};
      if (m_cfg.discRadialSegments > 1) {
        long double rStep = (discMaxR - discMinR) / m_cfg.discRadialSegments;
        radialValues.reserve(m_cfg.discRadialSegments);
        for (int ir = 0; ir <= m_cfg.discRadialSegments; ++ir) {
          radialValues.push_back(discMinR + ir * rStep);
        }
      } else {
        radialValues = {discMinR, discMaxR};
      }

      for (std::size_t ir = 1; ir < radialValues.size(); ++ir) {
        long double minR = radialValues[ir - 1];
        long double maxR = radialValues[ir];

        long double maxLocY = maxR * cosPhiHalf;
        long double minLocY = minR * cosPhiHalf;

        long double hR = 0.5 * (maxLocY + minLocY);
        long double hY = 0.5 * (maxLocY - minLocY);
        long double hXminY = minR * sinPhiHalf;
        long double hXmaxY = maxR * sinPhiHalf;

        auto tgTrapezoid =
            std::make_shared<Acts::TrapezoidBounds>(hXminY, hXmaxY, hY);

        for (int im = 0; im < m_cfg.discPhiSegments; ++im) {
          // Get the moduleTransform
          long double phi = -std::numbers::pi + im * phiStep;
          auto tgTransform = Transform3(
              Translation3(hR * std::cos(phi), hR * std::sin(phi),
                           discCenter.z()) *
              AngleAxis3(phi - std::numbers::pi / 2., Vector3::UnitZ()));

          // Create a new detector element per split
          auto tgDetectorElement = std::make_shared<Acts::TGeoDetectorElement>(
              tgIdentifier, tgNode, tgTransform, tgTrapezoid, tgThickness);

          tgDetectorElements.push_back(tgDetectorElement);
        }
      }

      return tgDetectorElements;
    }
  }

  // Cylinder segments are detected, attempt a split
  if (m_cfg.cylinderPhiSegments > 0 || m_cfg.cylinderLongitudinalSegments > 0) {
    if (sf.type() == Acts::Surface::Cylinder &&
        sf.bounds().type() == Acts::SurfaceBounds::eCylinder) {
      ACTS_DEBUG("- splitting detected for a Cylinder shaped sensor.");

      std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
          tgDetectorElements = {};
      tgDetectorElements.reserve(std::abs(m_cfg.cylinderPhiSegments) *
                                 std::abs(m_cfg.cylinderLongitudinalSegments));

      const auto& boundValues = sf.bounds().values();
      long double cylinderR = boundValues[Acts::CylinderBounds::eR];
      long double cylinderHalfZ =
          boundValues[Acts::CylinderBounds::eHalfLengthZ];

      long double phiStep = 2 * std::numbers::pi / m_cfg.cylinderPhiSegments;
      long double cosPhiHalf = std::cos(0.5 * phiStep);
      long double sinPhiHalf = std::sin(0.5 * phiStep);

      std::vector<long double> zValues = {};
      if (m_cfg.cylinderLongitudinalSegments > 1) {
        long double zStep =
            2 * cylinderHalfZ / m_cfg.cylinderLongitudinalSegments;
        zValues.reserve(m_cfg.cylinderLongitudinalSegments);
        for (int ir = 0; ir <= m_cfg.cylinderLongitudinalSegments; ++ir) {
          zValues.push_back(-cylinderHalfZ + ir * zStep);
        }
      } else {
        zValues = {-cylinderHalfZ, cylinderHalfZ};
      }

      long double planeR = cylinderR * cosPhiHalf;
      long double planeHalfX = cylinderR * sinPhiHalf;

      for (std::size_t iz = 1; iz < zValues.size(); ++iz) {
        long double minZ = zValues[iz - 1];
        long double maxZ = zValues[iz];

        long double planeZ = 0.5 * (minZ + maxZ);
        long double planeHalfY = 0.5 * (maxZ - minZ);

        auto tgRectangle =
            std::make_shared<Acts::RectangleBounds>(planeHalfX, planeHalfY);

        for (int im = 0; im < m_cfg.cylinderPhiSegments; ++im) {
          // Get the moduleTransform
          long double phi = -std::numbers::pi + im * phiStep;
          long double cosPhi = std::cos(phi);
          long double sinPhi = std::sin(phi);
          long double planeX = planeR * cosPhi;
          long double planeY = planeR * sinPhi;

          Acts::Vector3 planeCenter(planeX, planeY, planeZ);
          Acts::Vector3 planeAxisZ = {cosPhi, sinPhi, 0.};
          Acts::Vector3 planeAxisY{0., 0., 1.};
          Acts::Vector3 planeAxisX = planeAxisY.cross(planeAxisZ);

          RotationMatrix3 planeRotation;
          planeRotation.col(0) = planeAxisX;
          planeRotation.col(1) = planeAxisY;
          planeRotation.col(2) = planeAxisZ;

          // curvilinear surfaces are boundless
          Transform3 planeTransform{planeRotation};
          planeTransform.pretranslate(planeCenter);

          Transform3 tgTransform = planeTransform;

          // Create a new detector element per split
          auto tgDetectorElement = std::make_shared<Acts::TGeoDetectorElement>(
              tgIdentifier, tgNode, tgTransform, tgRectangle, tgThickness);

          tgDetectorElements.push_back(tgDetectorElement);
        }
      }
      return tgDetectorElements;
    }
  }
  return {tgde};
}
