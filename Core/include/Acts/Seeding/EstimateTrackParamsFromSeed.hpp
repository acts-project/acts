// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <vector>

namespace Acts {
/// @todo:
/// 1) Implement the simple Line and Circle fit based on Taubin Circle fit
/// 2) Implement the simple Line and Parabola fit (from HPS reconstruction by
/// Robert Johnson)

/// Estimate the track parameters on the xy plane from at least three space
/// points. It assumes the trajectory projection on the xy plane is a circle,
/// i.e. the magnetic field is along global z-axis.
///
/// The method is based on V. Karimaki NIM A305 (1991) 187-191:
/// https://doi.org/10.1016/0168-9002(91)90533-V
/// - no weights are used in Karimaki's fit, d0 is the distance of the point of
/// closest approach to the origin, 1/R is the curvature, phi is the angle of
/// the direction propagation (counter clockwise as positive) at the point of
/// cloest approach.
///
/// @tparam spacepoint_iterator_t The type of space point iterator
///
/// @param spBegin is the begin iterator for the space points
/// @param spEnd is the end iterator for the space points
///
/// @return optional estimated parameters d0, 1/R and phi

template <typename spacepoint_iterator_t>
std::optional<std::array<ActsScalar, 3>> estimateTrackParamsFromSeed(
    spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("estimateTrackParamsFromSeed", Logging::INFO));

  ActsScalar x2m = 0., xm = 0.;
  ActsScalar xym = 0.;
  ActsScalar y2m = 0., ym = 0.;
  ActsScalar r2m = 0., r4m = 0.;
  ActsScalar xr2m = 0., yr2m = 0.;

  size_t numSP = 0;
  for (spacepoint_iterator_t it = spBegin; it != spEnd; it++) {
    if (*it == nullptr) {
      continue;
    }

    const auto& sp = *it;

    ActsScalar x = sp->x();
    ActsScalar y = sp->y();
    ActsScalar r2 = x * x + y * y;
    x2m += x * x;
    xm += x;
    xym += x * y;
    y2m += y * y;
    ym += y;
    r2m += r2;
    r4m += r2 * r2;
    xr2m += x * r2;
    yr2m += y * r2;
    numSP++;
  }

  if (numSP < 3) {
    ACTS_ERROR("At least three valid space points are required.")
    return std::nullopt;
  }

  x2m = x2m / numSP;
  xm = xm / numSP;
  xym = xym / numSP;
  y2m = y2m / numSP;
  ym = ym / numSP;
  r2m = r2m / numSP;
  r4m = r4m / numSP;
  xr2m = xr2m / numSP;
  yr2m = yr2m / numSP;

  ActsScalar Cxx = x2m - xm * xm;
  ActsScalar Cxy = xym - xm * ym;
  ActsScalar Cyy = y2m - ym * ym;
  ActsScalar Cxr2 = xr2m - xm * r2m;
  ActsScalar Cyr2 = yr2m - ym * r2m;
  ActsScalar Cr2r2 = r4m - r2m * r2m;

  ActsScalar q1 = Cr2r2 * Cxy - Cxr2 * Cyr2;
  ActsScalar q2 = Cr2r2 * (Cxx - Cyy) - Cxr2 * Cxr2 + Cyr2 * Cyr2;

  ActsScalar phi = 0.5 * std::atan(2 * q1 / q2);
  ActsScalar k = (std::sin(phi) * Cxr2 - std::cos(phi) * Cyr2) * (1. / Cr2r2);
  ActsScalar delta = -k * r2m + std::sin(phi) * xm - std::cos(phi) * ym;

  ActsScalar rho = (2 * k) / (std::sqrt(1 - 4 * delta * k));
  ActsScalar d0 = (2 * delta) / (1 + std::sqrt(1 - 4 * delta * k));

  return std::array<ActsScalar, 3>{d0, rho, phi};
}

/// Estimate the full track parameters from three space points
///
/// This method is based on the conformal map transformation. It estimates the
/// full bound track parameters, i.e. (loc0, loc1, phi, theta, q/p, t) at the
/// bottom space point. The surface with the bottom space point on it should be
/// provided for the representation of the bound track parameters. The magnetic
/// field (which might be along any direction) is also necessary for the
/// momentum estimation.
///
/// It resembles the method used in ATLAS for the track parameters
/// estimated from seed, i.e. the function InDet::SiTrackMaker_xk::getAtaPlane
/// here:
/// https://acode-browser.usatlas.bnl.gov/lxr/source/athena/InnerDetector/InDetRecTools/SiTrackMakerTool_xk/src/SiTrackMaker_xk.cxx
///
/// @tparam spacepoint_iterator_t  The type of space point iterator
///
/// @param gctx is the geometry context
/// @param spBegin is the begin iterator for the space points
/// @param spEnd is the end iterator for the space points
/// @param surface is the surface of the bottom space point. The estimated bound
/// track parameters will be represented also at this surface
/// @param bField is the magnetic field vector
/// @param bFieldMin is the minimum magnetic field required to trigger the
/// estimation of q/pt
/// @param mass is the estimated particle mass
///
/// @return optional bound parameters
template <typename spacepoint_iterator_t>
std::optional<BoundVector> estimateTrackParamsFromSeed(
    const GeometryContext& gctx, spacepoint_iterator_t spBegin,
    spacepoint_iterator_t spEnd, const Surface& surface, Vector3 bField,
    ActsScalar bFieldMin = 0.1 * UnitConstants::T,
    ActsScalar mass = 139.57018 * UnitConstants::MeV) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("estimateTrackParamsFromSeed", Logging::INFO));

  // Convert bField to Tesla
  ActsScalar bFieldInTesla = bField.norm() / UnitConstants::T;
  ActsScalar bFieldMinInTesla = bFieldMin / UnitConstants::T;
  // Check if magnetic field is too small
  if (bFieldInTesla < bFieldMinInTesla) {
    // @todo shall we use straight-line estimation and use default q/pt in such
    // case?
    ACTS_WARNING("The magnetic field at the bottom space point: B = "
                 << bFieldInTesla << " T is smaller than |B|_min = "
                 << bFieldMinInTesla << " T. Estimation is not performed.")
    return std::nullopt;
  }

  // The global position of the bottom, middle and top space point
  Vector3 bGlobal, mGlobal, tGlobal;
  size_t numSP = 0;
  for (spacepoint_iterator_t it = spBegin; it != spEnd; it++) {
    if (*it == nullptr) {
      continue;
    }
    const auto& sp = *it;
    switch (numSP) {
      case 0: {
        bGlobal = Vector3(sp->x(), sp->y(), sp->z());
        break;
      }
      case 1: {
        mGlobal = Vector3(sp->x(), sp->y(), sp->z());
        break;
      }
      case 2: {
        tGlobal = Vector3(sp->x(), sp->y(), sp->z());
        break;
      }
      default: {
        ACTS_WARNING(
            "More than three space points are provided. The extra space point "
            "is not used for the "
            "estimation.");
      }
    }
    numSP++;
  }

  // Initialize the bound parameters vector
  BoundVector params = BoundVector::Zero();
  // Define a new coordinate frame with its origin at the bottom space point, z
  // axis long the magnetic field direction and y axis perpendicular to vector
  // from the bottom to middle space point. Hence, the projection of the middle
  // space point on the tranverse plane will be located at the x axis of the new
  // frame.
  Vector3 mRelVec = mGlobal - bGlobal;
  Vector3 newZAxis = bField.normalized();
  Vector3 newYAxis = newZAxis.cross(mRelVec).normalized();
  Vector3 newXAxis = newYAxis.cross(newZAxis);
  RotationMatrix3 rotation;
  rotation.col(0) = newXAxis;
  rotation.col(1) = newYAxis;
  rotation.col(2) = newZAxis;
  // The center of the new frame is at the bottom space point
  Translation3 trans(bGlobal);
  // The transform which constructs the new frame
  Transform3 transform(trans * rotation);

  // The coordinate of the middle and top space point in the new frame
  Vector3 local1 = transform.inverse() * mGlobal;
  Vector3 local2 = transform.inverse() * tGlobal;

  // Lambda to transform the coordinates to the (u, v) space
  auto uvTransform = [](const Vector3& local) -> Vector2 {
    Vector2 uv;
    ActsScalar denominator = local.x() * local.x() + local.y() * local.y();
    uv.x() = local.x() / denominator;
    uv.y() = local.y() / denominator;
    return uv;
  };
  // The uv1.y() should be zero
  Vector2 uv1 = uvTransform(local1);
  Vector2 uv2 = uvTransform(local2);

  // A,B are slope and intercept of the straight line in the u,v plane
  // connecting the three points.
  ActsScalar A = uv2.y() / (uv2.x() - uv1.x());
  ActsScalar B = uv2.y() - A * uv2.x();
  // Curvature (with a sign) estimate
  ActsScalar rho = -2.0 * B / std::hypot(1., A);
  // The projection of the top space point on the transverse plane of the new
  // frame
  ActsScalar rn = local2.x() * local2.x() + local2.y() * local2.y();
  // The (1/tanTheta) of momentum in the new frame, corrected for curvature
  // effects (@note why 0.04?)
  ActsScalar invTanTheta =
      local2.z() * std::sqrt(1. / rn) / (1. + 0.04 * rho * rho * rn);
  // The momentum direction in the new frame (the center of the circle has the
  // coordinate (-1.*A/(2*B), 1./(2*B)))
  Vector3 transDirection(1., A, std::hypot(1, A) * invTanTheta);
  // Transform it back to the original frame
  Vector3 direction = rotation * transDirection.normalized();

  // The phi and theta
  params[eBoundPhi] = VectorHelpers::phi(direction);
  params[eBoundTheta] = VectorHelpers::theta(direction);

  // Transform the bottom space point to local coordinates
  auto lpResult = surface.globalToLocal(gctx, bGlobal, direction);
  if (not lpResult.ok()) {
    ACTS_ERROR("Global to local transformation did not succeed.");
    return std::nullopt;
  }
  Vector2 bLocal = lpResult.value();
  // The estimated loc0 and loc1
  params[eBoundLoc0] = bLocal.x();
  params[eBoundLoc1] = bLocal.y();

  // The estimated q/pt in [GeV/c]^-1 (note that the pt is the projection of
  // momentum in the transverse plane of the new frame)
  ActsScalar qOverPt = rho * (UnitConstants::m) / (0.3 * bFieldInTesla);
  // The estimated q/p in [GeV/c]^-1
  params[eBoundQOverP] = qOverPt / std::hypot(1., invTanTheta);

  // The estimated momentum and projection of the momentum in the magnetic field
  // diretion
  ActsScalar pInGeV = std::abs(1.0 / params[eBoundQOverP]);
  ActsScalar pzInGeV = 1.0 / std::abs(qOverPt) * invTanTheta;
  ActsScalar massInGeV = mass / UnitConstants::GeV;
  // The velocity along the magnetic field diretion (i.e. z axis of the new
  // frame)
  ActsScalar vz = pzInGeV / std::hypot(pInGeV, massInGeV);
  // The z coordinate of the bottom space point along the magnetic field
  // direction
  ActsScalar pathz = bGlobal.dot(bField) / bField.norm();
  // The estimated time
  params[eBoundTime] = pathz / vz;

  return params;
}

}  // namespace Acts
