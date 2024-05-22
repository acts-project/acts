
// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelSurfaceConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoTrap.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoMaterial.h>
#include <GeoModelKernel/GeoShape.h>

Acts::GeoModelSensitiveSurface
Acts::GeoModelSurfaceConverter::convertToSensitiveSurface(
    const GeoFullPhysVol& geoPhysVol) {
  const GeoLogVol* logVol = geoPhysVol.getLogVol();

  const Transform3& transform = geoPhysVol.getAbsoluteTransform(nullptr);
  auto fullPhysVol = dynamic_cast<const GeoFullPhysVol*>(&geoPhysVol);

  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();
    if (geoShape != nullptr) {
      Vector3 translation = transform.translation();
      RotationMatrix3 rotation = transform.rotation();

      Transform3 surfaceTransform = Acts::Transform3::Identity();
      surfaceTransform.translation() = translation;

      // Try if its a box shape
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

        // Create the surface
        ActsScalar halfX = halfLengths[xIndex];
        ActsScalar halfY = halfLengths[yIndex];
        // Create the surface
        auto rectangleBounds =
            std::make_shared<Acts::RectangleBounds>(halfX, halfY);

        auto detectorElement = std::make_shared<Acts::GeoModelDetectorElement>(
            geoPhysVol, rectangleBounds, surfaceTransform,
            2 * halfLengths[zIndex]);
        auto surface = detectorElement->surface().getSharedPtr();
        return std::make_tuple(detectorElement, surface);
      }

      // Try if it is trapezoid of type GeoTrap
      auto geoTrap = dynamic_cast<const GeoTrap*>(geoShape);
      if (geoTrap != nullptr) {
        // Create the surface
        ActsScalar minHalfX = geoTrap->getDxdyndzn();
        ActsScalar maxHalfX = geoTrap->getDxdypdzn();
        ActsScalar halfY    = geoTrap->getDydzp();
        auto trapezoidBounds =
            std::make_shared<Acts::TrapezoidBounds>(minHalfX, maxHalfX, halfY);

        auto detectorElement = std::make_shared<Acts::GeoModelDetectorElement>(
            geoPhysVol, trapezoidBounds, surfaceTransform,
            2 * geoTrap->getZHalfLength());
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
        ActsScalar diffX = std::abs(halfX2-halfX1);
        ActsScalar diffY = std::abs(halfY2-halfY1);

        // And its interpretation
        ActsScalar minHalfX = diffX > diffY ? halfX1 : halfY1;
        ActsScalar maxHalfX = diffX > diffY ? halfX2 : halfY2;
        ActsScalar thickness = diffX > diffY ? diffY : diffX;

        // Adjust the rotation matrix 
        RotationMatrix3 trotation = surfaceTransform.rotation();
        if (diffX > diffY){
          // Rotation is x, z, y ... acyclic, hence the sign change
          trotation.col(1) = surfaceTransform.rotation().col(2);
          trotation.col(2) = -surfaceTransform.rotation().col(1);
        } else {
          // Rotation is y, z, x ... cyclic, hence no sign change
          trotation.col(0) = surfaceTransform.rotation().col(1);
          trotation.col(1) = surfaceTransform.rotation().col(2);
          trotation.col(2) = surfaceTransform.rotation().col(0);
          surfaceTransform.linear() = trotation;
        }

        auto trapezoidBounds = std::make_shared<Acts::TrapezoidBounds>(
             minHalfX, maxHalfX, halfZ);

        auto detectorElement = std::make_shared<Acts::GeoModelDetectorElement>(
            geoPhysVol, trapezoidBounds, surfaceTransform, thickness);
        auto surface = detectorElement->surface().getSharedPtr();
        return std::make_tuple(detectorElement, surface);
      }
    }
  }

  return std::make_tuple(nullptr, nullptr);
}
