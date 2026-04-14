// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/GeoModel/detail/GeoTubeConverter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/GeoModel/GeoModelConversionError.hpp"

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/Units.h>

using namespace Acts;

Result<ActsPlugins::GeoModelSensitiveSurface>
ActsPlugins::detail::GeoTubeConverter::operator()(
    const PVConstLink& geoPV, const GeoTube& geoTube,
    const Transform3& absTransform, SurfaceBoundFactory& boundFactory,
    bool sensitive) const {
  /// auto-calculate the unit length conversion
  static constexpr double unitLength =
      UnitConstants::mm / GeoModelKernelUnits::millimeter;

  // Create the surface transform
  Transform3 transform = Transform3::Identity();
  transform.translation() = unitLength * absTransform.translation();
  transform.linear() = absTransform.rotation();

  // Create the surface
  double innerRadius = unitLength * geoTube.getRMin();
  double outerRadius = unitLength * geoTube.getRMax();
  double halfZ = unitLength * geoTube.getZHalfLength();

  if (targetShape == Surface::SurfaceType::Straw) {
    // Create the element and the surface
    auto lineBounds = boundFactory.makeBounds<LineBounds>(outerRadius, halfZ);
    if (!sensitive) {
      auto surface = Surface::makeShared<StrawSurface>(transform, lineBounds);
      return std::make_tuple(nullptr, surface);
    }

    auto surf0 = Surface::makeShared<StrawSurface>(transform, lineBounds);
    auto detEl0 = GeoModelDetectorElement::createDetectorElement(
        geoPV, transform, 2 * outerRadius, surf0);
    return std::make_tuple(detEl0, surf0);
    // Next option is translation to disc
  } else if (targetShape == Surface::SurfaceType::Disc) {
    auto radialBounds =
        std::make_shared<RadialBounds>(innerRadius, outerRadius);
    if (!sensitive) {
      auto surface = Surface::makeShared<DiscSurface>(transform, radialBounds);
      return std::make_tuple(nullptr, surface);
    }

    // Create the element and the surface
    auto surf1 = Surface::makeShared<DiscSurface>(transform, radialBounds);
    auto detEl1 = GeoModelDetectorElement::createDetectorElement(
        geoPV, transform, 2 * halfZ, surf1);
    return std::make_tuple(detEl1, surf1);
  }
  // Finally cylinder to cylinder
  auto cylinderBounds = std::make_shared<CylinderBounds>(outerRadius, halfZ);
  if (!sensitive) {
    auto surface =
        Surface::makeShared<CylinderSurface>(transform, cylinderBounds);
    return std::make_tuple(nullptr, surface);
  }
  // Create the element and the surface
  auto surf2 = Surface::makeShared<CylinderSurface>(transform, cylinderBounds);
  auto detEl2 = GeoModelDetectorElement::createDetectorElement(
      geoPV, transform, outerRadius - innerRadius, surf2);
  return std::make_tuple(detEl2, surf2);
}
