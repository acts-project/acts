// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>

namespace Acts {

namespace {

Transform3 estimationFrameLocalToGlobal(const Vector3& sp0, const Vector3& sp1,
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
  const Translation3 translation(sp0);
  // The transform which constructs the new frame
  return translation * rotation;
}

double computeDzDs(double A, double B, const Vector3& local0,
                   const Vector3& local2) {
  const auto computeLocalPhi = [&](const Vector2& local) -> double {
    // Scaled radius vector from circle center
    const Vector2 r = 2 * B * local - Vector2(-A, 1);

    return std::atan2(r.y(), r.x());
  };

  const double localPhi0 = computeLocalPhi(local0.head<2>());
  const double localPhi2 = computeLocalPhi(local2.head<2>());

  const double dZ = local2.z() - local0.z();
  const double dPhi = localPhi2 - localPhi0;

  // Apply the sinc correction to account for the fact that the particle do not
  // follow a straight line in the R-Z plane. This is especially important for
  // high delta phi and/or strongly bent tracks.
  const double sincCorrection = sinc(dPhi / 2);

  const double dzds =
      sincCorrection * dZ / (local2.head<2>() - local0.head<2>()).norm();

  return dzds;
}

struct ConformalMappingResult {
  Vector2 uv1;
  Vector2 uv2;
  Vector2 duv;
  double A;
  double B;
  double bOverS;
  double dzds;
};

ConformalMappingResult performConformalMapping(const Vector3& local1,
                                               const Vector3& local2) {
  ConformalMappingResult r{};
  r.uv1 = local1.head<2>() / local1.head<2>().squaredNorm();
  r.uv2 = local2.head<2>() / local2.head<2>().squaredNorm();
  r.duv = r.uv2 - r.uv1;
  r.A = r.duv.y() / r.duv.x();
  r.B = r.uv1.y() - r.A * r.uv1.x();
  r.bOverS = (r.uv1.y() * r.uv2.x() - r.uv2.y() * r.uv1.x()) / r.duv.norm();
  r.dzds = computeDzDs(r.A, r.B, local1, local2);
  return r;
}

Vector3 computeLocalTangent(const ConformalMappingResult& cm,
                            const Vector2& local) {
  // Scaled radius vector from circle center
  const Vector2 r = 2 * cm.B * local - Vector2(-cm.A, 1);

  // Tangent perpendicular to radius
  const Vector3 t(-r.y(), r.x(), r.norm() * cm.dzds);

  return t.normalized();
}

}  // namespace

}  // namespace Acts

Acts::FreeVector Acts::estimateTrackParamsFromSeed(const Vector3& sp0,
                                                   const Vector3& sp1,
                                                   const Vector3& sp2,
                                                   const Vector3& bField) {
  return estimateTrackParamsFromSeed(sp0, 0, sp1, sp2, bField);
}

Acts::FreeVector Acts::estimateTrackParamsFromSeed(
    const Vector3& sp0, const double t0, const Vector3& sp1, const Vector3& sp2,
    const Vector3& bField, Vector3* tangent0, Vector3* tangent1,
    Vector3* tangent2) {
  const Transform3 transform = estimationFrameLocalToGlobal(sp0, sp1, bField);

  // Local coordinates
  const Vector3 local0 = Vector3::Zero();
  const Vector3 local1 = transform.inverse() * sp1;
  const Vector3 local2 = transform.inverse() * sp2;

  // Conformal mapping
  const ConformalMappingResult cm = performConformalMapping(local1, local2);

  // The tangent vector at the bottom space point
  const Vector3 direction =
      transform.linear() * computeLocalTangent(cm, local0.head<2>());

  // Initialize the free parameters vector
  FreeVector params = FreeVector::Zero();

  // The bottom space point position
  params.segment<3>(eFreePos0) = sp0;

  // The estimated direction
  params.segment<3>(eFreeDir0) = direction;

  // The estimated q/pt in [GeV/c]^-1 (note that the pt is the projection of
  // momentum on the transverse plane of the new frame)
  const double qOverPt = 2 * cm.bOverS / bField.norm();
  // The estimated q/p in [GeV/c]^-1
  params[eFreeQOverP] = qOverPt / fastHypot(1, cm.dzds);

  // The time parameter is set to the time of the bottom space point
  params[eFreeTime] = t0;

  if (tangent0 != nullptr) {
    *tangent0 = direction;
  }
  if (tangent1 != nullptr) {
    *tangent1 = transform.linear() * computeLocalTangent(cm, local1.head<2>());
  }
  if (tangent2 != nullptr) {
    *tangent2 = transform.linear() * computeLocalTangent(cm, local2.head<2>());
  }

  return params;
}

Acts::Result<Acts::BoundVector> Acts::estimateTrackParamsFromSeed(
    const GeometryContext& gctx, const Surface& surface, const Vector3& sp0,
    const double t0, const Vector3& sp1, const Vector3& sp2,
    const Vector3& bField) {
  const FreeVector freeParams =
      estimateTrackParamsFromSeed(sp0, t0, sp1, sp2, bField);
  return transformFreeToBoundParameters(freeParams, surface, gctx);
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
