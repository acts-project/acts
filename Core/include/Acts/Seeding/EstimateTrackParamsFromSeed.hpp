// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
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
#include <iostream>
#include <optional>
#include <vector>

namespace Acts {
/// @todo:
/// 1) Implement the simple Line and Circle fit based on Taubin Circle fit
/// 2) Implement the simple Line and Parabola fit (from HPS reconstruction by
/// Robert
/// Johnson)

/// Estimate the track parameters on transverse plane
///
/// @note Using method based on V. Karimaki NIM A305 (1991) 187-191 - no weights
/// are used in Karimaki's fit, d0 is the distance of the closest approach to
/// the origin, 1/R is the curvature, phi is the angle of the direction
/// propagation (counter clockwise as positive)
///
/// @tparam external_spacepoint_t The type of space point
/// @param sps is a vector of space points with positions in unit of mm
///
/// @return optional estimated parameters d0, 1/R and phi
template <typename external_spacepoint_t>
std::optional<std::array<double, 3>> estimateTrackParamsFromSeed(
    const std::vector<external_spacepoint_t*>& sps) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("estimateTrackParamsFromSeed", Logging::INFO));
  if (sps.empty()) {
    ACTS_ERROR("No space points exsit.")
    return std::nullopt;
  }

  double x2m = 0., xm = 0.;
  double xym = 0.;
  double y2m = 0., ym = 0.;
  double r2m = 0., r4m = 0.;
  double xr2m = 0., yr2m = 0.;

  size_t numSP = sps.size();

  for (const auto& sp : sps) {
    double x = sp->x();
    double y = sp->y();
    double r2 = x * x + y * y;
    x2m += x * x;
    xm += x;
    xym += x * y;
    y2m += y * y;
    ym += y;
    r2m += r2;
    r4m += r2 * r2;
    xr2m += x * r2;
    yr2m += y * r2;
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

  double Cxx = x2m - xm * xm;
  double Cxy = xym - xm * ym;
  double Cyy = y2m - ym * ym;
  double Cxr2 = xr2m - xm * r2m;
  double Cyr2 = yr2m - ym * r2m;
  double Cr2r2 = r4m - r2m * r2m;

  double q1 = Cr2r2 * Cxy - Cxr2 * Cyr2;
  double q2 = Cr2r2 * (Cxx - Cyy) - Cxr2 * Cxr2 + Cyr2 * Cyr2;

  double phi = 0.5 * std::atan(2 * q1 / q2);
  double k = (std::sin(phi) * Cxr2 - std::cos(phi) * Cyr2) * (1. / Cr2r2);
  double delta = -k * r2m + std::sin(phi) * xm - std::cos(phi) * ym;

  double rho = (2 * k) / (std::sqrt(1 - 4 * delta * k));
  double d = (2 * delta) / (1 + std::sqrt(1 - 4 * delta * k));

  return std::array<double, 3>{d, rho, phi};
}

///  Estimate the full track parameters from space points
///
/// @note This resembles the method used in ATLAS for the track parameters
/// estimated from seed, i.e. the function InDet::SiTrackMaker_xk::getAtaPlane
/// here:
/// https://acode-browser.usatlas.bnl.gov/lxr/source/athena/InnerDetector/InDetRecTools/SiTrackMakerTool_xk/src/SiTrackMaker_xk.cxx
///
/// @note This method gives an estimate of the track parameters of a seed using
/// a conformal map transformation The track parameters are of the form loc1,
/// loc2, phi, theta, q/p. phi0 is the angle of the track direction with respect
/// the origin, positive when counter clock-wise
///
/// @tparam external_spacepoint_t The type of space point
///
/// @param gctx is the geometry context
/// @param sps is the vector of space points with positions in unit of mm
/// @param surface is the surface at which the bound track parameters will be
/// represented
/// @param bFieldZ is the magnetic field in z direction in T
/// @param bFieldZMin is the minimum magnetic field in z direction in T when
/// estimation of q/p is done
/// @param mass is the estimated particle mass
///
/// @return optional bound parameters
template <typename external_spacepoint_t>
std::optional<BoundVector> estimateTrackParamsFromSeed(
    const GeometryContext& gctx, const std::vector<external_spacepoint_t>& sps,
    const Surface& surface, double bFieldZ,
    double bFieldZMin = 0.1 * UnitConstants::T,
    double mass = 139.57018 * UnitConstants::MeV) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("estimateTrackParamsFromSeed", Logging::INFO));
  if (sps.size() < 3) {
    ACTS_ERROR("At least three space points are required.")
    return std::nullopt;
  }

  // Make sure the pt and bFieldZ are in expected units
  double bFieldZInTesla = bFieldZ / UnitConstants::T;
  double bFieldZMinInTesla = bFieldZMin / UnitConstants::T;
  // Check if magnetic field is too small
  if (bFieldZInTesla < bFieldZMinInTesla) {
    // @note shall we use straight-line estimation and use default q/pt in such
    // case?
    ACTS_WARNING("The magnetic field at the first space point: Bz = "
               << bFieldZInTesla << " T is too small.")
    return std::nullopt;
  }

  // Initialize the bound parameters vector
  BoundVector params = BoundVector::Zero();

  double x0 = sps[0]->x();
  double y0 = sps[0]->y();
  double z0 = sps[0]->z();
  // Define the locations of the second and third space points with respect to
  // the first one in the Space point vector
  double x1 = sps[1]->x() - x0;
  double y1 = sps[1]->y() - y0;
  double x2 = sps[2]->x() - x0;
  double y2 = sps[2]->y() - y0;
  double z2 = sps[2]->z() - z0;

  // Define conformal map variables
  double u1 = 1. / std::sqrt(x1 * x1 + y1 * y1);
  // denominator for conformal mapping
  double rn = x2 * x2 + y2 * y2;
  double r2 = 1. / rn;
  // coordinate system for conformal mapping
  double a = x1 * u1;
  double b = y1 * u1;
  // (u, v) coordinates of third SP in conformal mapping
  double u2 = (a * x2 + b * y2) * r2;
  double v2 = (a * y2 - b * x2) * r2;
  // A,B are slope and intercept of the straight line in the u,v plane
  // connecting the three points
  double A = v2 / (u2 - u1);
  double B = 2. * (v2 - A * u2);
  // Curvature estimate : (2R)²=(1+A²)/b² => 1/2R = b/sqrt(1+A²) = B /
  // sqrt(1+A²)
  // @todo checks the sign of rho
  double rho = -1.0 * B / std::sqrt(1. + A * A);

  // Estimate of the track dz/dr (1/tanTheta), corrected for curvature effects
  // @note why 0.04?
  double invTanTheta = z2 * std::sqrt(r2) / (1. + 0.04 * rho * rho * rn);
  // The estimated phi
  params[2] = std::atan2(b + a * A, a - b * A);
  // The estimated theta
  params[3] = std::atan2(1., invTanTheta);
  // The estimated q/pt in [GeV/c]^-1 assuming the space
  // point global positions have units in mm
  double estQOverPt = rho * 1000 / (0.3 * bFieldZInTesla);
  // The estimated q/p in [GeV/c]^-1
  params[4] = estQOverPt / std::sqrt(1. + invTanTheta * invTanTheta);
  // The estimated pz and p in GeV/c
  double pzInGeV = 1.0 / std::abs(estQOverPt) * invTanTheta;
  double pInGeV = std::abs(1.0 / params[4]);
  double massInGeV = mass / UnitConstants::GeV;
  double vz = pzInGeV / std::sqrt(pInGeV * pInGeV + massInGeV * massInGeV);
  // The estimated time in ms?
  params[5] = sps[0]->z() / vz;

  // The estimated momentum direction
  double sinTheta = std::sin(params[3]);
  double cosTheta = std::cos(params[3]);
  double sinPhi = std::sin(params[2]);
  double cosPhi = std::cos(params[2]);
  Vector3 direction(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
  // Transform the first space point to the provided surface
  Vector3 global(x0, y0, z0);
  auto lpResult = surface.globalToLocal(gctx, global, direction);
  if (not lpResult.ok()) {
    ACTS_ERROR("Global to local transformation did not succeed.");
    return std::nullopt;
  }
  Vector2 local = lpResult.value();
  // The estimated loc0 and loc1
  params[0] = local[0];
  params[1] = local[1];

  return params;
}

}  // namespace Acts
