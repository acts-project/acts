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
#include "Acts/Seeding/detail/CircleFit.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Eigenvalues>

namespace Acts {

namespace {

Transform3 estimationFrameLocalToGlobal(const Vector3& sp0, const Vector3& sp1,
                                        const Vector3& bField) {
  // Define a new coordinate frame with its origin at the bottom space point, z
  // axis long the magnetic field direction and y axis perpendicular to vector
  // from the bottom to middle space point. Hence, the projection of the middle
  // space point on the transverse plane will be located at the x axis of the
  // new frame.
  //
  // The construction is robust against a vanishing field (the z axis falls back
  // to the global z axis) and against a reference vector parallel to the field
  // (any orthogonal y axis will do, as the transverse orientation does not
  // affect the global result).
  const Vector3 relVec = sp1 - sp0;
  const double bMag = bField.norm();
  const Vector3 newZAxis = (bMag > std::numeric_limits<double>::epsilon())
                               ? Vector3(bField / bMag)
                               : Vector3(Vector3::UnitZ());
  Vector3 newYAxis = newZAxis.cross(relVec);
  if (newYAxis.norm() < std::numeric_limits<double>::epsilon()) {
    newYAxis = newZAxis.unitOrthogonal();
  }
  newYAxis.normalize();
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

namespace {

class TrackParamsEstimationErrorCategory : public std::error_category {
 public:
  const char* name() const noexcept final {
    return "TrackParamsEstimationError";
  }

  std::string message(int c) const final {
    using Acts::TrackParamsEstimationError;

    switch (static_cast<TrackParamsEstimationError>(c)) {
      case TrackParamsEstimationError::NotEnoughSpacePoints:
        return "At least three space points are required";
      case TrackParamsEstimationError::DegenerateFit:
        return "The space point configuration leads to a degenerate fit";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(Acts::TrackParamsEstimationError e) {
  static TrackParamsEstimationErrorCategory c;
  return {static_cast<int>(e), c};
}

Acts::Result<Acts::FreeVector> Acts::estimateTrackParamsFromSpacePoints(
    std::span<const Vector3> spacePoints, const Vector3& bField, double t0,
    std::size_t geometricRefineIterations, std::span<const double> weights) {
  if (spacePoints.size() < 3) {
    return Result<FreeVector>::failure(
        TrackParamsEstimationError::NotEnoughSpacePoints);
  }
  assert((weights.empty() || weights.size() == spacePoints.size()) &&
         "weights must be empty or match the space points");

  // Per-point weight accessor: an empty span means uniform weights of one.
  const auto w = [&](std::size_t i) {
    return weights.empty() ? 1. : weights[i];
  };

  const Vector3& reference = spacePoints.front();

  // Estimation frame: field along the local +z axis, fixed by the first two
  // points and the field. All points are expressed in this frame.
  const Transform3 transform =
      estimationFrameLocalToGlobal(reference, spacePoints[1], bField);
  const Transform3 toLocal = transform.inverse();
  std::vector<Vector3> local;
  local.reserve(spacePoints.size());
  for (const Vector3& sp : spacePoints) {
    local.push_back(toLocal * sp);
  }

  FreeVector params = FreeVector::Zero();
  params.segment<3>(eFreePos0) = reference;
  params[eFreeTime] = t0;

  const double bMag = bField.norm();

  // Algebraic circle fit in the transverse plane, optionally refined
  // geometrically. Both stages use the per-point weights.
  detail::CircleFit circle = detail::fitCircleTaubin(local, weights);
  if (circle.valid && geometricRefineIterations > 0) {
    detail::refineCircleGeometric(circle, local, geometricRefineIterations,
                                  weights);
  }

  if (!circle.valid) {
    // Straight-line limit (collinear transverse points, i.e. a straight track
    // or a vanishing field): the direction is the principal axis of the local
    // points and the momentum is dropped (q/p stays zero).
    double sumW = 0.;
    Vector3 mean = Vector3::Zero();
    for (std::size_t i = 0; i < local.size(); ++i) {
      sumW += w(i);
      mean += w(i) * local[i];
    }
    if (!(sumW > 0.)) {
      return Result<FreeVector>::failure(
          TrackParamsEstimationError::DegenerateFit);
    }
    mean /= sumW;
    SquareMatrix3 cov = SquareMatrix3::Zero();
    for (std::size_t i = 0; i < local.size(); ++i) {
      const Vector3 d = local[i] - mean;
      cov += w(i) * d * d.transpose();
    }
    Eigen::SelfAdjointEigenSolver<SquareMatrix3> solver(cov);
    // Eigenvalues are ascending; the largest measures the spread along the
    // principal axis. A non-positive value means the points coincide and no
    // direction can be recovered.
    if (solver.info() != Eigen::Success || (solver.eigenvalues()[2] <= 0)) {
      return Result<FreeVector>::failure(
          TrackParamsEstimationError::DegenerateFit);
    }
    Vector3 localDir = solver.eigenvectors().col(2);
    if (localDir.dot(local.back() - local.front()) < 0.) {
      localDir = -localDir;
    }
    localDir.normalize();
    params.segment<3>(eFreeDir0) = transform.linear() * localDir;
    return Result<FreeVector>::success(params);
  }

  const Vector2 center = circle.center;
  const double radius = circle.radius;

  const Vector2 refXy = local.front().head<2>();
  const Vector2 secondXy = local[1].head<2>();

  // Tangent at the reference point: the +90 degree rotation of the radial
  // vector (from the center to the reference), oriented along the travel
  // direction (reference -> second point).
  const Vector2 radial = refXy - center;
  Vector2 tangent(-radial.y(), radial.x());
  const double rotSense = (tangent.dot(secondXy - refXy) >= 0.) ? 1. : -1.;
  tangent = (rotSense * tangent).normalized();

  // Charge sign: the Lorentz force q v x B points towards the circle center for
  // the true charge sign. With B along the local +z axis, v x B = (v_y, -v_x).
  const double toCenterProj = tangent.y() * (center.x() - refXy.x()) -
                              tangent.x() * (center.y() - refXy.y());
  const double qSign = (toCenterProj >= 0.) ? 1. : -1.;

  // Exact R-Z line fit: z is linear in the transverse arc length
  // s = radius * (turning angle about the center), measured along travel. The
  // turning angle is unwrapped point-to-point (consecutive points are assumed
  // less than half a turn apart).
  const double phiRef =
      std::atan2(refXy.y() - center.y(), refXy.x() - center.x());
  double phiPrev = phiRef;
  double sumW = 0., sumS = 0., sumZ = 0., sumSS = 0., sumSZ = 0.;
  for (std::size_t i = 0; i < local.size(); ++i) {
    const Vector3& p = local[i];
    double phi = std::atan2(p.y() - center.y(), p.x() - center.x());
    while (phi - phiPrev > std::numbers::pi) {
      phi -= 2. * std::numbers::pi;
    }
    while (phi - phiPrev < -std::numbers::pi) {
      phi += 2. * std::numbers::pi;
    }
    phiPrev = phi;
    const double wi = w(i);
    const double s = rotSense * radius * (phi - phiRef);
    const double z = p.z();
    sumW += wi;
    sumS += wi * s;
    sumZ += wi * z;
    sumSS += wi * s * s;
    sumSZ += wi * s * z;
  }
  const double denom = sumW * sumSS - sumS * sumS;
  const double lambda =
      (std::abs(denom) > std::numeric_limits<double>::epsilon())
          ? (sumW * sumSZ - sumS * sumZ) / denom
          : 0.;

  // Assemble the direction: the local tangent is (t_x, t_y, lambda) because s
  // is the transverse arc length, then rotate back into the global frame.
  Vector3 localDir(tangent.x(), tangent.y(), lambda);
  localDir.normalize();
  params.segment<3>(eFreeDir0) = transform.linear() * localDir;

  // Momentum: p_T[GeV] = |B|_native * R[mm], hence q/p_T = qSign / (R * |B|).
  // Without a field the momentum cannot be estimated, so q/p is left at zero.
  if (bMag > std::numeric_limits<double>::epsilon()) {
    const double qOverPt = qSign / (radius * bMag);
    params[eFreeQOverP] = qOverPt / fastHypot(1., lambda);
  }

  return Result<FreeVector>::success(params);
}

Acts::Result<Acts::BoundVector> Acts::estimateTrackParamsFromSpacePoints(
    const GeometryContext& gctx, const Surface& surface,
    std::span<const Vector3> spacePoints, const Vector3& bField, double t0,
    std::size_t geometricRefineIterations, std::span<const double> weights) {
  Result<FreeVector> freeParams = estimateTrackParamsFromSpacePoints(
      spacePoints, bField, t0, geometricRefineIterations, weights);
  if (!freeParams.ok()) {
    return Result<BoundVector>::failure(freeParams.error());
  }
  return transformFreeToBoundParameters(*freeParams, surface, gctx);
}
