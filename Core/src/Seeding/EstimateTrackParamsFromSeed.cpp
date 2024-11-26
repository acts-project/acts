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
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/VoidNavigator.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <numbers>

namespace {

class EstimateTrackParamsFromSeedErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category.
  const char* name() const noexcept final {
    return "EstimateTrackParamsFromSeedError";
  }

  // Return what each enum means in text.
  std::string message(int c) const final {
    using Acts::EstimateTrackParamsFromSeedError;

    switch (static_cast<EstimateTrackParamsFromSeedError>(c)) {
      case EstimateTrackParamsFromSeedError::SurfaceIntersectionFailed:
        return "Surface intersection failed";
      case EstimateTrackParamsFromSeedError::PropagationFailed:
        return "Propagation failed";
      default:
        return "unknown";
    }
  }
};

}  // namespace

std::error_code Acts::make_error_code(
    Acts::EstimateTrackParamsFromSeedError e) {
  static EstimateTrackParamsFromSeedErrorCategory c;
  return {static_cast<int>(e), c};
}

Acts::FreeVector Acts::estimateTrackParamsFromSeed(const Vector3& sp0,
                                                   const Vector3& sp1,
                                                   const Vector3& sp2,
                                                   const Vector3& bField) {
  // Define a new coordinate frame with its origin at the bottom space point, z
  // axis long the magnetic field direction and y axis perpendicular to vector
  // from the bottom to middle space point. Hence, the projection of the middle
  // space point on the transverse plane will be located at the x axis of the
  // new frame.
  Vector3 relVec = sp1 - sp0;
  Vector3 newZAxis = bField.normalized();
  Vector3 newYAxis = newZAxis.cross(relVec).normalized();
  Vector3 newXAxis = newYAxis.cross(newZAxis);
  RotationMatrix3 rotation;
  rotation.col(0) = newXAxis;
  rotation.col(1) = newYAxis;
  rotation.col(2) = newZAxis;
  // The center of the new frame is at the bottom space point
  Translation3 trans(sp0);
  // The transform which constructs the new frame
  Transform3 transform(trans * rotation);

  // The coordinate of the middle and top space point in the new frame
  Vector3 local1 = transform.inverse() * sp1;
  Vector3 local2 = transform.inverse() * sp2;

  // In the new frame the bottom sp is at the origin, while the middle
  // sp in along the x axis. As such, the x-coordinate of the circle is
  // at: x-middle / 2.
  // The y coordinate can be found by using the straight line passing
  // between the mid point between the middle and top sp and perpendicular to
  // the line connecting them
  Vector2 circleCenter;
  circleCenter(0) = 0.5 * local1(0);

  ActsScalar deltaX21 = local2(0) - local1(0);
  ActsScalar sumX21 = local2(0) + local1(0);
  // straight line connecting the two points
  // y = a * x + c (we don't care about c right now)
  // we simply need the slope
  // we compute 1./a since this is what we need for the following computation
  ActsScalar ia = deltaX21 / local2(1);
  // Perpendicular line is then y = -1/a *x + b
  // we can evaluate b given we know a already by imposing
  // the line passes through P = (0.5 * (x2 + x1), 0.5 * y2)
  ActsScalar b = 0.5 * (local2(1) + ia * sumX21);
  circleCenter(1) = -ia * circleCenter(0) + b;
  // Radius is a signed distance between circleCenter and first sp, which is at
  // (0, 0) in the new frame. Sign depends on the slope a (positive vs negative)
  int sign = ia > 0 ? -1 : 1;
  const ActsScalar R = circleCenter.norm();
  ActsScalar invTanTheta =
      local2.z() / (2 * R * std::asin(local2.head<2>().norm() / (2 * R)));
  // The momentum direction in the new frame (the center of the circle has the
  // coordinate (-1.*A/(2*B), 1./(2*B)))
  ActsScalar A = -circleCenter(0) / circleCenter(1);
  Vector3 transDirection(1., A, fastHypot(1, A) * invTanTheta);
  // Transform it back to the original frame
  Vector3 direction = rotation * transDirection.normalized();

  // Initialize the free parameters vector
  FreeVector params = FreeVector::Zero();

  // The bottom space point position
  params.segment<3>(eFreePos0) = sp0;

  // The estimated direction
  params.segment<3>(eFreeDir0) = direction;

  // The estimated q/pt in [GeV/c]^-1 (note that the pt is the projection of
  // momentum on the transverse plane of the new frame)
  ActsScalar qOverPt = sign / (bField.norm() * R);
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

      // transverse momentum contribution
      variance += std::pow(config.initialSigmaPtRel * params[eBoundQOverP], 2);

      // theta contribution
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

Acts::Result<Acts::BoundVector> Acts::estimateTrackParamsFromSeedAtSurface(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const Surface& surface, const Vector3& sp0, const Vector3& sp1,
    const Vector3& sp2,
    const std::shared_ptr<const MagneticFieldProvider>& bField,
    const BasePropagator* propagator, const Acts::Logger& logger) {
  MagneticFieldProvider::Cache bFieldCache = bField->makeCache(mctx);

  Result<Vector3> bFieldResult = bField->getField(sp0, bFieldCache);
  if (!bFieldResult.ok()) {
    ACTS_INFO("The magnetic field failed.");
    return Result<BoundVector>::failure(bFieldResult.error());
  }
  const Vector3& bFieldVector = bFieldResult.value();

  FreeVector freeParams =
      estimateTrackParamsFromSeed(sp0, sp1, sp2, bFieldVector);

  Vector3 origin = sp0;
  Vector3 direction = freeParams.segment<3>(eFreeDir0);
  double qOverPt = freeParams[eFreeQOverP];

  auto lpResult = surface.globalToLocal(gctx, origin, direction);
  if (!lpResult.ok()) {
    // no cov transport matrix is needed here
    // particle hypothesis does not matter here
    CurvilinearTrackParameters estimatedParams(freeParams.segment<4>(eFreePos0),
                                               direction, qOverPt, std::nullopt,
                                               ParticleHypothesis::pion());

    auto surfaceIntersection =
        surface.intersect(gctx, origin, direction).closest();

    if (!surfaceIntersection.isValid()) {
      ACTS_INFO(
          "The surface does not intersect with the origin and estimated "
          "direction.");
      // TODO
      return Result<BoundVector>::failure(
          EstimateTrackParamsFromSeedError::SurfaceIntersectionFailed);
    }

    Direction propagatorDirection =
        Direction::fromScalarZeroAsPositive(surfaceIntersection.pathLength());

    PropagatorPlainOptions propagatorOptions(gctx, mctx);
    propagatorOptions.direction = propagatorDirection;

    Result<BoundTrackParameters> result =
        (propagator != nullptr)
            ? propagator->propagateToSurface(estimatedParams, surface,
                                             propagatorOptions)
            : Propagator(EigenStepper<>(bField), VoidNavigator(),
                         logger().cloneWithSuffix("Propagator"))
                  .propagateToSurface(estimatedParams, surface,
                                      propagatorOptions);
    if (!result.ok()) {
      ACTS_INFO("The propagation failed.");
      return Result<BoundVector>::failure(
          EstimateTrackParamsFromSeedError::PropagationFailed);
    }

    return Result<BoundVector>::success(result.value().parameters());
  }

  BoundVector params = BoundVector::Zero();
  params[eBoundPhi] = VectorHelpers::phi(direction);
  params[eBoundTheta] = VectorHelpers::theta(direction);
  params[eBoundQOverP] = freeParams[eFreeQOverP];

  Vector2 bottomLocalPos = lpResult.value();
  // The estimated loc0 and loc1
  params[eBoundLoc0] = bottomLocalPos.x();
  params[eBoundLoc1] = bottomLocalPos.y();
  // We need the bottom space point time later
  params[eBoundTime] = 0;

  return Result<BoundVector>::success(params);
}
