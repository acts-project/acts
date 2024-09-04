// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/detail/GeoPolygonConverter.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConversionError.hpp"
#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/Units.h>

Acts::Result<Acts::GeoModelSensitiveSurface>
Acts::detail::GeoPolygonConverter::operator()(
    const PVConstLink& geoPV, const GeoSimplePolygonBrep& polygon,
    const Transform3& absTransform, bool sensitive) const {
  /// auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / GeoModelKernelUnits::millimeter;

  // Create the surface transform
  Transform3 transform = Transform3::Identity();
  transform.translation() = unitLength * absTransform.translation();
  auto rotation = absTransform.rotation();
  // Get the half lengths
  int nVertices = polygon.getNVertices();
  std::vector<std::vector<double>> vertices;

  for (int i = 0; i < nVertices; i++) {
    vertices.push_back({polygon.getXVertex(i), polygon.getYVertex(i)});
  }
  // sort based on the y-coordinate
  std::sort(vertices.begin(), vertices.end(),
            [](const std::vector<double>& a, const std::vector<double>& b) {
              return a[1] < b[1];
            });
  if (nVertices == 4) {
    double hlxnegy = fabs(vertices[0][0] - vertices[1][0]) / 2;
    double hlxposy = fabs(vertices[2][0] - vertices[3][0]) / 2;
    double hly = fabs(vertices[0][1] - vertices[3][1]) / 2;
    std::vector<ActsScalar> halfLengths = {hlxnegy, hlxposy, hly};

    // Create the surface
    Vector3 colX = rotation.col(0);
    Vector3 colY = rotation.col(1);
    Vector3 colZ = rotation.col(2);
    rotation.col(0) = colX;
    rotation.col(1) = colY;
    rotation.col(2) = colZ;
    transform.linear() = rotation;
    // Create the surface bounds
    ActsScalar halfXnegY = unitLength * halfLengths[0];
    ActsScalar halfXposY = unitLength * halfLengths[1];
    ActsScalar halfY = unitLength * halfLengths[2];
    auto trapBounds =
        std::make_shared<Acts::TrapezoidBounds>(halfXnegY, halfXposY, halfY);
    if (!sensitive) {
      auto surface = Surface::makeShared<PlaneSurface>(transform, trapBounds);
      return std::make_tuple(nullptr, surface);
    }
    // Create the element and the surface
    auto detectorElement =
        GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
            geoPV, trapBounds, transform, unitLength * polygon.getDZ());
    auto surface = detectorElement->surface().getSharedPtr();
    // Return the detector element and surface
    return std::make_tuple(detectorElement, surface);
  } else if (nVertices == 6) {
    double hlxnegy = fabs(vertices[0][0] - vertices[1][0]) / 2;
    double hlxzeroy = fabs(vertices[2][0] - vertices[3][0]) / 2;
    double hlxposy = fabs(vertices[4][0] - vertices[5][0]) / 2;
    double hly = fabs(vertices[0][1] - vertices[4][1]) / 2;
    std::vector<ActsScalar> halfLengths = {hlxnegy, hlxzeroy, hlxposy, hly,
                                           hly};

    // Create the surface

    Vector3 colX = rotation.col(0);
    Vector3 colY = rotation.col(1);
    Vector3 colZ = rotation.col(2);
    rotation.col(0) = colX;
    rotation.col(1) = colY;
    rotation.col(2) = colZ;
    transform.linear() = rotation;

    // Create the surface bounds
    ActsScalar halfXnegY = unitLength * halfLengths[0];
    ActsScalar halfXzeroY = unitLength * halfLengths[1];
    ActsScalar halfXposY = unitLength * halfLengths[2];
    ActsScalar halfYnegX = unitLength * halfLengths[3];
    ActsScalar halfYposX = unitLength * halfLengths[4];

    auto diamondBounds = std::make_shared<Acts::DiamondBounds>(
        halfXnegY, halfXzeroY, halfXposY, halfYnegX, halfYposX);
    if (!sensitive) {
      auto surface =
          Surface::makeShared<PlaneSurface>(transform, diamondBounds);
      return std::make_tuple(nullptr, surface);
    }
    // Create the element and the surface
    auto detectorElement =
        GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
            geoPV, diamondBounds, transform, unitLength * polygon.getDZ());
    auto surface = detectorElement->surface().getSharedPtr();
    // Return the detector element and surface
    return std::make_tuple(detectorElement, surface);
  } else {
    throw std::runtime_error("GeoSimplePolygonBrep with " +
                             std::to_string(nVertices) +
                             " can not be converted");
  }
}
