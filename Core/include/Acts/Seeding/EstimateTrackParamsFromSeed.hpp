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
#include "Acts/Seeding/Seed.hpp"
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
  if (sps.empty()) {
    ACTS_LOCAL_LOGGER(
        getDefaultLogger("TrackParametersEstimation", Logging::INFO));
    ACTS_FATAL("No space points exsit.")
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
/// @param sps is the vector of space points with positions in unit of mm
/// @param transform is the transform which helps to represent the track
/// parameters
/// @param bFieldZ is the magnetic field in z direction in T
/// @param ptMin is the minimum transverse momentum allowed in GeV
/// @param ptTolerance is the tolerance of pt in percent
/// @param mass is the estimated particle mass
///
/// @return optional bound parameters
template <typename external_spacepoint_t>
std::optional<BoundVector> estimateTrackParamsFromSeed(
    const std::vector<external_spacepoint_t>& sps, const Transform3& transform,
    double bFieldZ, double ptMin = 0.5 * UnitConstants::GeV,
    double ptTolerance = 0.1, double mass = 139.57018 * UnitConstants::MeV) {
  if (sps.size() < 3) {
    ACTS_LOCAL_LOGGER(
        getDefaultLogger("TrackParametersEstimation", Logging::INFO));
    ACTS_FATAL("At least three space points are required.")
    return std::nullopt;
  }

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

  // Project to the surface
  const auto matrix = transform.matrix();
  Vector3 Ax = matrix.block<3, 1>(0, 0);
  Vector3 Ay = matrix.block<3, 1>(0, 1);
  // Center of the surface in the global frame
  Vector3 D = matrix.block<3, 1>(0, 3);
  // Location of the first SP w.r.t. the surface center
  Vector3 d(x0 - D[0], y0 - D[1], z0 - D[2]);

  // Initialize the bound parameters vector
  BoundVector params = BoundVector::Zero();
  // The estimated loc0 and loc1
  params[0] = d.dot(Ax);
  params[1] = d.dot(Ay);

  // Make sure the ptMin and bFieldZ are in expected units
  double ptMinInGeV = ptMin / UnitConstants::GeV;
  double bFieldZInTesla = bFieldZ / UnitConstants::T;

  // The 1/tanTheta inverse and qOverPt
  double invTanTheta;
  double qOverPt;
  // Check if magnetic field is more than 0.1 Tesla
  if (bFieldZInTesla > 0.1) {
    // Estimate of the track dz/dr (1/tanTheta), corrected for curvature effects
    invTanTheta = z2 * std::sqrt(r2) / (1. + 0.04 * rho * rho * rn);
    // The estimated phi
    params[2] = std::atan2(b + a * A, a - b * A);
    // The estimated inverse transverse momentum in [GeV/c]^-1 assumes the space
    // point global positions have units in mm
    // @note the curvature needs to be in the unit of [m]^-1.
    qOverPt = rho * 1000 / (0.3 * bFieldZInTesla);
  } else {
    // Use straight-line estimate (no curvature correction)
    invTanTheta = z2 * sqrt(r2);
    // The estimated phi
    params[2] = std::atan2(y2, x2);
    // No pt estimation, assume min pt and positive charge
    qOverPt = 1. / ptMinInGeV;
  }

  // Return nullopt if the estimate pt is not not within requirement
  if (std::abs(qOverPt) * ptMinInGeV > (1. + ptTolerance)) {
    return std::nullopt;
  }

  // The estimated theta
  params[3] = std::atan2(1., invTanTheta);
  // The estimated q/p in [GeV/c]^-1 using q/pt and 1.0/tanTheta
  params[4] = qOverPt / std::sqrt(1. + invTanTheta * invTanTheta);

  // The estimated pz and p in GeV/c
  double pzInGeV = 1.0 / std::abs(qOverPt) * invTanTheta;
  double pInGeV = 1.0 / params[4];
  double massInGeV = mass / UnitConstants::GeV;
  double vz = pzInGeV / std::sqrt(pInGeV * pInGeV + massInGeV * massInGeV);
  // The estimated time in ms?
  params[5] = sps[0]->z() / vz;
  return params;
}

}  // namespace Acts
