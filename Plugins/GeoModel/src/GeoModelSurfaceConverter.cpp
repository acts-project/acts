
// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelSurfaceConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <limits>

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/GeoTrap.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/Units.h>

Acts::GeoModelSensitiveSurface
Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(
    const GeoFullPhysVol& geoPhysVol, const Options& options) {
  // Get the logical volume first
  const GeoLogVol* logVol = geoPhysVol.getLogVol();
  const Transform3& transform = geoPhysVol.getAbsoluteTransform(nullptr);
  /// auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / GeoModelKernelUnits::millimeter;

  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();
    if (geoShape != nullptr) {
      Vector3 translation = unitLength * transform.translation();
      RotationMatrix3 rotation = transform.rotation();

      Transform3 surfaceTransform = Acts::Transform3::Identity();
      surfaceTransform.translation() = translation;

      // (1) Try if its a box shape
      auto geoBox = dynamic_cast<const GeoBox*>(geoShape);
      if (geoBox != nullptr) {
        std::vector<ActsScalar> halfLengths = {geoBox->getXHalfLength(),
                                               geoBox->getYHalfLength(),
                                               geoBox->getZHalfLength()};
        // Create the surface
        auto minElement =
            std::min_element(halfLengths.begin(), halfLengths.end());
        auto zIndex = std::distance(halfLengths.begin(), minElement);
        std::size_t yIndex = zIndex > 0u ? zIndex - 1u : 2u;
        std::size_t xIndex = yIndex > 0u ? yIndex - 1u : 2u;

        Vector3 colX = rotation.col(xIndex);
        Vector3 colY = rotation.col(yIndex);
        Vector3 colZ = rotation.col(zIndex);
        rotation.col(0) = colX;
        rotation.col(1) = colY;
        rotation.col(2) = colZ;
        surfaceTransform.linear() = rotation;

        // Create the surface bounds
        ActsScalar halfX = unitLength * halfLengths[xIndex];
        ActsScalar halfY = unitLength * halfLengths[yIndex];
        auto rectangleBounds =
            std::make_shared<Acts::RectangleBounds>(halfX, halfY);
        // Create the element and the surface
        auto detectorElement =
            GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
                geoPhysVol, rectangleBounds, surfaceTransform,
                2 * unitLength * halfLengths[zIndex]);

        auto surface = detectorElement->surface().getSharedPtr();
        return std::make_tuple(detectorElement, surface);
      }

      // (2) Try if it is trapezoid of type GeoTrap
      auto geoTrap = dynamic_cast<const GeoTrap*>(geoShape);
      if (geoTrap != nullptr) {
        // Create the traepzoidal bounds
        ActsScalar minHalfX = unitLength * geoTrap->getDxdyndzn();
        ActsScalar maxHalfX = unitLength * geoTrap->getDxdypdzn();
        ActsScalar halfY = unitLength * geoTrap->getDydzp();
        auto trapezoidBounds =
            std::make_shared<Acts::TrapezoidBounds>(minHalfX, maxHalfX, halfY);
        // Create the element and the surface
        auto detectorElement =
            GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
                geoPhysVol, trapezoidBounds, surfaceTransform,
                2 * unitLength * geoTrap->getZHalfLength());
        auto surface = detectorElement->surface().getSharedPtr();
        return std::make_tuple(detectorElement, surface);
      }

      // Try if it is a trapezoid of type GeoTrd
      auto geoTrd = dynamic_cast<const GeoTrd*>(geoShape);
      if (geoTrd != nullptr) {
        // GeoTrd is defined that halfZ needs to map onto surface halfY
        // Create the surface
        ActsScalar halfX1 = geoTrd->getXHalfLength1();
        ActsScalar halfX2 = geoTrd->getXHalfLength2();
        ActsScalar halfY1 = geoTrd->getYHalfLength1();
        ActsScalar halfY2 = geoTrd->getYHalfLength2();
        ActsScalar halfZ = geoTrd->getZHalfLength();

        // The diffs
        ActsScalar diffX = std::abs(halfX2 - halfX1);
        ActsScalar diffY = std::abs(halfY2 - halfY1);

        // If both X and Y are trapezoidal - consider
        if (diffX > 2 * std::numeric_limits<ActsScalar>::epsilon() &&
            diffY > 2 * std::numeric_limits<ActsScalar>::epsilon()) {
          throw std::invalid_argument(
              "GeoModelSurfaceConverter: GeoTrd conversion to Trapezoid "
              "ambiguous.");
        }

        // And its interpretation
        ActsScalar minHalfX = unitLength * (diffX > diffY ? halfX1 : halfY1);
        ActsScalar maxHalfX = unitLength * (diffX > diffY ? halfX2 : halfY2);
        ActsScalar thickness = unitLength * (diffX > diffY ? diffY : diffX);

        // This is a convetion of the TrapezoidBounds
        int swapZ = (maxHalfX < minHalfX) ? -1 : 1;
        if (swapZ < 0) {
          std::swap(minHalfX, maxHalfX);
        }

        // Adjust the rotation matrix
        RotationMatrix3 trotation = surfaceTransform.rotation();
        if (diffX > diffY) {
          // Rotation is x, z, y ... acyclic, hence the sign change
          trotation.col(1) = swapZ * surfaceTransform.rotation().col(2);
          trotation.col(2) = -swapZ * surfaceTransform.rotation().col(1);
        } else {
          // Rotation is y, z, x ... cyclic, hence no sign change
          trotation.col(0) = surfaceTransform.rotation().col(1);
          trotation.col(1) = swapZ * surfaceTransform.rotation().col(2);
          trotation.col(2) = swapZ * surfaceTransform.rotation().col(0);
        }
        surfaceTransform.linear() = trotation;
        auto trapezoidBounds =
            std::make_shared<Acts::TrapezoidBounds>(minHalfX, maxHalfX, halfZ);
        // Create the element and the surface
        auto detectorElement =
            GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
                geoPhysVol, trapezoidBounds, surfaceTransform, thickness);
        auto surface = detectorElement->surface().getSharedPtr();
        return std::make_tuple(detectorElement, surface);
      }

      // (3) Try if it is a tube
      auto geoTube = dynamic_cast<const GeoTube*>(geoShape);
      if (geoTube != nullptr) {
        // Create the surface
        ActsScalar innerRadius = unitLength * geoTube->getRMin();
        ActsScalar outerRadius = unitLength * geoTube->getRMax();
        ActsScalar halfZ = unitLength * geoTube->getZHalfLength();

        if (options.cylinderTargetType == Acts::Surface::SurfaceType::Straw) {
          // Create the element and the surface
          auto lineBounds =
              std::make_shared<Acts::LineBounds>(outerRadius, halfZ);
          auto detectorElement =
              GeoModelDetectorElement::createDetectorElement<StrawSurface>(
                  geoPhysVol, lineBounds, surfaceTransform, 2 * outerRadius);
          auto surface = detectorElement->surface().getSharedPtr();
          return std::make_tuple(detectorElement, surface);
          // Next option is tranlsation to disc
        } else if (options.cylinderTargetType ==
                   Acts::Surface::SurfaceType::Disc) {
          auto radialBounds =
              std::make_shared<Acts::RadialBounds>(innerRadius, outerRadius);
          // Create the element and the surface
          auto detectorElement =
              GeoModelDetectorElement::createDetectorElement<DiscSurface>(
                  geoPhysVol, radialBounds, surfaceTransform, 2 * halfZ);
          auto surface = detectorElement->surface().getSharedPtr();
          return std::make_tuple(detectorElement, surface);
        }
        // Finally cylinder to cylinder
        auto cylinderBounds =
            std::make_shared<Acts::CylinderBounds>(outerRadius, halfZ);
        // Create the element and the surface
        auto detectorElement =
            GeoModelDetectorElement::createDetectorElement<CylinderSurface>(
                geoPhysVol, cylinderBounds, surfaceTransform,
                outerRadius - innerRadius);
        auto surface = detectorElement->surface().getSharedPtr();
        return std::make_tuple(detectorElement, surface);
      }
    }
  }

  return std::make_tuple(nullptr, nullptr);
}
