// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

Acts::FreeVector Acts::estimateTrackParamsFromSeed(const Vector3& sp0,
                                                   const Vector3& sp1,
                                                   const Vector3& sp2,
                                                   const Vector3& bField) {
  // Define a new coordinate frame with its origin at the bottom space point, z
  // axis long the magnetic field direction and y axis perpendicular to vector
  // from the bottom to middle space point. Hence, the projection of the middle
  // space point on the transverse plane will be located at the x axis of the
  // new frame.
  const Vector3 relVec = sp1 - sp0;
  const Vector3 newZAxis = bField.normalized();
  const Vector3 newYAxis = newZAxis.cross(relVec).normalized();
  const Vector3 newXAxis = newYAxis.cross(newZAxis);
  RotationMatrix3 rotation;
  rotation.col(0) = newXAxis;
  rotation.col(1) = newYAxis;
  rotation.col(2) = newZAxis;
  // The center of the new frame is at the bottom space point
  const Translation3 trans(sp0);
  // The transform which constructs the new frame
  const Transform3 transform(trans * rotation);

  // The coordinate of the middle and top space point in the new frame
  const Vector3 local1 = transform.inverse() * sp1;
  const Vector3 local2 = transform.inverse() * sp2;

  // We determine the center of the circle directly given the 3 points and using
  // Cramer's rule.
  Vector2 circleCenter = Vector2::Zero();
  {
    const double det = 2 * (local1.x() * local2.y() - local2.x() * local1.y());

    // Check if points are aligned and return straight line estimate if so
    if (std::abs(det) < 1e-12) {
      FreeVector params = FreeVector::Zero();
      params.segment<3>(eFreePos0) = sp0;
      params.segment<3>(eFreeDir0) = (sp2 - sp0).normalized();
      params[eFreeQOverP] = 0;
      return params;
    }

    const double z1 = local1.head<2>().squaredNorm();
    const double z2 = local2.head<2>().squaredNorm();

    const double nom1 = local1.y() * z2 - local2.y() * z1;
    const double nom2 = local2.x() * z1 - local1.x() * z2;

    circleCenter.x() = nom1 / det;
    circleCenter.y() = nom2 / det;
  }
  const int sign = circleCenter.y() >= 0 ? 1 : -1;
  const double R = circleCenter.norm();

  const double invTanTheta =
      local2.z() / (2 * R * std::asin(local2.head<2>().norm() / (2 * R)));
  // The momentum direction in the new frame (the center of the circle has the
  // coordinate (-1.*A/(2*B), 1./(2*B)))
  const double A = -circleCenter.x() / circleCenter.y();
  const Vector3 transDirection(1., A, fastHypot(1, A) * invTanTheta);
  // Transform it back to the original frame
  const Vector3 direction = rotation * transDirection.normalized();

  // Initialize the free parameters vector
  FreeVector params = FreeVector::Zero();

  // The bottom space point position
  params.segment<3>(eFreePos0) = sp0;

  // The estimated direction
  params.segment<3>(eFreeDir0) = direction;

  // The estimated q/pt in [GeV/c]^-1 (note that the pt is the projection of
  // momentum on the transverse plane of the new frame)
  const double qOverPt = sign / (bField.norm() * R);
  // The estimated q/p in [GeV/c]^-1
  params[eFreeQOverP] = qOverPt / fastHypot(1., invTanTheta);

  return params;
}

Acts::BoundMatrix Acts::estimateTrackParamCovariance(
    const EstimateTrackParamCovarianceConfig& config, const BoundVector& params,
    bool hasTime) {
  assert((params[eBoundTheta] > 0 && params[eBoundTheta] < std::numbers::pi) &&
         "Theta must be in the range (0, pi)");

  BoundSquareMatrix result = BoundSquareMatrix::Zero();

  for (std::size_t i = eBoundLoc0; i < eBoundSize; ++i) {
    double sigma = config.initialSigmas[i];
    double variance = sigma * sigma;

    if (i == eBoundQOverP) {
      // note that we rely on the fact that sigma theta is already computed
      double varianceTheta = result(eBoundTheta, eBoundTheta);

      // contribution from sigma(q/pt)
      variance += std::pow(
          config.initialSigmaQoverPt * std::sin(params[eBoundTheta]), 2);

      // contribution from sigma(pt)/pt
      variance += std::pow(config.initialSigmaPtRel * params[eBoundQOverP], 2);

      // contribution from sigma(theta)
      variance +=
          varianceTheta *
          std::pow(params[eBoundQOverP] / std::tan(params[eBoundTheta]), 2);
    }

    if (i == eBoundTime && !hasTime) {
      // Inflate the time uncertainty if no time measurement is available
      variance *= config.noTimeVarInflation;
    }

    // Inflate the initial covariance
    variance *= config.initialVarInflation[i];

    result(i, i) = variance;
  }

  return result;
}
