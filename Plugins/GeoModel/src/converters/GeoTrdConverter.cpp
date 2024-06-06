// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/converters/GeoTrdConverter.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConversionError.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/Units.h>

Acts::Result<Acts::GeoModelSensitiveSurface>
Acts::detail::GeoTrdConverter::operator()(const GeoFullPhysVol& geoFPV,
                                          const GeoTrd& geoTrd,
                                          const Transform3& absTransform,
                                          bool sensitive) const {
  /// auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / GeoModelKernelUnits::millimeter;

  // Create the surface transform
  Transform3 transform = Transform3::Identity();
  transform.translation() = unitLength * absTransform.translation();

  // GeoTrd is defined that halfZ needs to map onto surface halfY
  // Create the surface
  ActsScalar halfX1 = geoTrd.getXHalfLength1();
  ActsScalar halfX2 = geoTrd.getXHalfLength2();
  ActsScalar halfY1 = geoTrd.getYHalfLength1();
  ActsScalar halfY2 = geoTrd.getYHalfLength2();
  ActsScalar halfZ = geoTrd.getZHalfLength();

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

  // This is a convention of the TrapezoidBounds
  int swapZ = (maxHalfX < minHalfX) ? -1 : 1;
  if (swapZ < 0) {
    std::swap(minHalfX, maxHalfX);
  }

  // Adjust the rotation matrix
  RotationMatrix3 trotation = absTransform.rotation();

  if (diffX > diffY) {
    // Rotation is x, z, y ... acyclic, hence the sign change
    trotation.col(1) = swapZ * absTransform.rotation().col(2);
    trotation.col(2) = -swapZ * absTransform.rotation().col(1);
  } else {
    // Rotation is y, z, x ... cyclic, hence no sign change
    trotation.col(0) = absTransform.rotation().col(1);
    trotation.col(1) = swapZ * absTransform.rotation().col(2);
    trotation.col(2) = swapZ * absTransform.rotation().col(0);
  }
  transform.linear() = trotation;

  auto trapezoidBounds =
      std::make_shared<TrapezoidBounds>(minHalfX, maxHalfX, halfZ);
  // std::cout << "     TrapezoidBounds: minHalfX=" << minHalfX << ", maxHalfX="
  // << maxHalfX << ", halfz=" << halfZ << std::endl;
  if (!sensitive) {
    auto surface =
        Surface::makeShared<PlaneSurface>(transform, trapezoidBounds);
    return std::make_tuple(nullptr, surface);
  }

  // Create the element and the surface
  auto detectorElement =
      GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
          geoFPV, trapezoidBounds, transform, thickness);
  auto surface = detectorElement->surface().getSharedPtr();
  return std::make_tuple(detectorElement, surface);
}
