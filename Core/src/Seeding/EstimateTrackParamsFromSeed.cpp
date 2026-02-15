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
  // Define a new coordinate frame with its origin at the bottom spacepoint, z
  // axis long the magnetic field direction and y axis perpendicular to vector
  // from the bottom to middle spacepoint. Hence, the projection of the middle
  // spacepoint on the transverse plane will be located at the x axis of the
  // new frame.
  const Vector3 relVec = sp1 - sp0;
  const Vector3 newZAxis = bField.normalized();
  const Vector3 newYAxis = newZAxis.cross(relVec).normalized();
  const Vector3 newXAxis = newYAxis.cross(newZAxis);
  RotationMatrix3 rotation;
  rotation.col(0) = newXAxis;
  rotation.col(1) = newYAxis;
  rotation.col(2) = newZAxis;
  // The center of the new frame is at the bottom spacepoint
  const Translation3 trans(sp0);
  // The transform which constructs the new frame
  const Transform3 transform(trans * rotation);

  // The coordinate of the middle and top spacepoint in the new frame
  const Vector3 local1 = transform.inverse() * sp1;
  const Vector3 local2 = transform.inverse() * sp2;

  // Use the uv-plane to estimate the circle parameters
  const Vector2 uv1 = local1.head<2>() / local1.head<2>().squaredNorm();
  const Vector2 uv2 = local2.head<2>() / local2.head<2>().squaredNorm();
  const Vector2 deltaUV2 = uv2 - uv1;
  const double A = deltaUV2.y() / deltaUV2.x();
  const double bOverS =
      (uv1.y() * uv2.x() - uv2.y() * uv1.x()) / deltaUV2.norm();

  const double invTanTheta = local2.z() / local2.head<2>().norm();
  const Vector3 transDirection(1, A, fastHypot(1, A) * invTanTheta);

  // Transform it back to the original frame
  const Vector3 direction = rotation * transDirection.normalized();

  // Initialize the free parameters vector
  FreeVector params = FreeVector::Zero();

  // The bottom spacepoint position
  params.segment<3>(eFreePos0) = sp0;

  // The estimated direction
  params.segment<3>(eFreeDir0) = direction;

  // The estimated q/pt in [GeV/c]^-1 (note that the pt is the projection of
  // momentum on the transverse plane of the new frame)
  const double qOverPt = 2 * bOverS / bField.norm();
  // The estimated q/p in [GeV/c]^-1
  params[eFreeQOverP] = qOverPt / fastHypot(1, invTanTheta);

  return params;
}

Acts::BoundMatrix Acts::estimateTrackParamCovariance(
    const EstimateTrackParamCovarianceConfig& config, const BoundVector& params,
    bool hasTime) {
  assert((params[eBoundTheta] > 0 && params[eBoundTheta] < std::numbers::pi) &&
         "Theta must be in the range (0, pi)");

  BoundMatrix result = BoundMatrix::Zero();

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
