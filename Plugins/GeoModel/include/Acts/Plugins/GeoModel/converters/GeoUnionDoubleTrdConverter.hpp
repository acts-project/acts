// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoShiftConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Utilities/Result.hpp"
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/LineBounds.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>
#include <Acts/Surfaces/TrapezoidBounds.hpp>

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoShapeUnion.h>
#include <GeoModelKernel/GeoTrd.h>

#include "GeoTrdConverter.hpp"
#include "GeoShiftConverter.hpp"
#include "GeoBoxConverter.hpp"

namespace Acts {

namespace detail {

  auto distanceLinePoint(const Acts::Vector3 &lineA, const Acts::Vector3 &lineB, const Acts::Vector3 &p) {
    auto dir = lineB - lineA;
    auto ap = p - lineA;
    return ap.cross(dir).norm() / dir.norm();
  }

  auto distanceParallelLines(const Acts::Vector3 &pointA, const Acts::Vector3 &pointB, const Acts::Vector3 &dir) {
    return std::sqrt( (pointA - pointB).cross(dir).norm() / dir.norm() );
  }

  std::shared_ptr<Surface> unionIsBigTrapezoid(const Surface &a, const Surface &b, double gapTolerance = 0.2) {
    // Assume the following Layout for this:
    //
    //  0 ______________________ 1
    //    \                    /
    //     \        B         /
    //      \________________/
    //     3                 2
    //      0 _____________ 1
    //        \           /
    //         \    A    /
    //          \_______/
    //          3       2

    if( a.bounds().type() != SurfaceBounds::eTrapezoid || b.bounds().type() != SurfaceBounds::eTrapezoid ) {
      return nullptr;
    }

    const auto &boundsA = static_cast<const TrapezoidBounds &>(a.bounds());
    const auto &boundsB = static_cast<const TrapezoidBounds &>(b.bounds());

    auto pa = a.polyhedronRepresentation({}, 0);
    auto pb = b.polyhedronRepresentation({}, 0);

    const auto &vtxsa = pa.vertices;
    const auto &vtxsb = pb.vertices;

    auto printv = [&](const auto &vs) {
    for(const auto &v : vs) {
      std::cout << "[ " << v.transpose() << " ]  ";
    }
    };

    std::cout << "vtxs A: ";
    printv(vtxsa);
    std::cout << std::endl;

    std::cout << "vtxs B: ";
    printv(vtxsb);
    std::cout << std::endl;

    const auto d1 = (vtxsb[3] - vtxsa[0]).norm();
    const auto d2 = (vtxsb[2] - vtxsa[1]).norm();
    if( d1 > gapTolerance || d2 > gapTolerance ) {
      std::cout << "Gap tolerance violated" << std::endl;
    }
    std::cout << "Gap is small enough " << d1 << " " << d2 << std::endl;

    auto P1 = vtxsa[0] + 0.5 * (vtxsb[3] - vtxsa[0]);
    auto dist1 = distanceLinePoint(vtxsa[3], vtxsb[0], P1);
    std::cout << "distance point to line: " << dist1 << std::endl;

    auto P2 = vtxsa[1] + 0.5 * (vtxsb[2] - vtxsa[1]);
    auto dist2 = distanceLinePoint(vtxsa[2], vtxsb[1], P2);
    std::cout << "distance point to line: " << dist2 << std::endl;

    std::cout << "vector A32 = " << (vtxsa[2]-vtxsa[3]).normalized().transpose() << std::endl;
    std::cout << "vector A01 = " << (vtxsa[1]-vtxsa[0]).normalized().transpose() << std::endl;
    std::cout << "vector B01 = " << (vtxsb[1]-vtxsb[0]).normalized().transpose() << std::endl;
    std::cout << "vector B32 = " << (vtxsb[2]-vtxsb[3]).normalized().transpose() << std::endl;


    const Vector3 dir = vtxsa[1]-vtxsa[0];
    const auto distanceTrapezoids = distanceParallelLines(vtxsa[0], vtxsb[3], dir);

    std::cout << "distance trapezoids: " << distanceTrapezoids << std::endl;

    auto hlxny = boundsA.values()[TrapezoidBounds::eHalfLengthXposY];
    auto hlxpy = boundsB.values()[TrapezoidBounds::eHalfLengthXnegY];

    Acts::Vector3 mp1 = vtxsa[3] + 0.5 * (vtxsa[2] - vtxsa[3]);
    Acts::Vector3 mp2 = vtxsb[0] + 0.5 * (vtxsb[1] - vtxsb[0]);

    auto hly = 0.5 * (mp2 - mp1).norm();

    std::cout << "pars hlxny: " << hlxny << ", hlxpy: " << hlxpy << ", hly: " << hly << std::endl;

    auto trapezoidBounds =
      std::make_shared<TrapezoidBounds>(hlxpy, hlxny, hly);
    auto transform = a.transform({});
    transform.translate(Vector3{0.f, -(hly - boundsA.values()[TrapezoidBounds::eHalfLengthY]), 0.f});

    return Surface::makeShared<PlaneSurface>(transform, trapezoidBounds);
  }
}

/// This is super hacky and ugly
/// TODO Please to not merge!!!
struct GeoUnionConverter : public IGeoShapeConverter {
  bool useA = true;

  auto convert(const GeoFullPhysVol& geoFPV, const GeoShape& shape, const Transform3& transform,
               bool sensitive) const -> Result<GeoModelSensitiveSurface> {
    if (auto trd = dynamic_cast<const GeoTrd*>(&shape); trd != nullptr) {
      return detail::GeoTrdConverter{}(geoFPV, *trd, transform, sensitive);
    }
    if (auto box = dynamic_cast<const GeoBox*>(&shape); box != nullptr) {
      return detail::GeoBoxConverter{}(geoFPV, *box, transform, sensitive);
    }
    if (auto shift = dynamic_cast<const GeoShapeShift*>(&shape); shift != nullptr) {
      return detail::GeoShiftConverter{}(geoFPV, *shift, transform, sensitive);
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }

  Acts::Result<Acts::GeoModelSensitiveSurface> toSensitiveSurface(
      const GeoFullPhysVol& geoFPV) const override {
    // Retrieve logcal volume and absolute transform
    const GeoLogVol* logVol = geoFPV.getLogVol();
    const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
    if (logVol != nullptr) {
      const GeoShape* geoShape = logVol->getShape();

      // This should be only called when this is clear
      // As I said, super hacky
      auto concreteShape = static_cast<const GeoShapeUnion*>(geoShape);

      auto pick = useA ? concreteShape->getOpA() : concreteShape->getOpB();

      if (concreteShape != nullptr) {
        return convert(geoFPV, *pick, transform, true);
      }
      return Result<GeoModelSensitiveSurface>::failure(
          GeoModelConversionError::WrongShapeForConverter);
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::MissingLogicalVolume);
  }

  Acts::Result<std::shared_ptr<Acts::Surface>> toPassiveSurface(
      const GeoFullPhysVol&) const override {
    throw std::runtime_error("not implemented");
  }
};

}  // namespace Acts
