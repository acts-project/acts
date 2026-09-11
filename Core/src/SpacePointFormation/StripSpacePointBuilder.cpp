// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/SpacePointFormation/PixelSpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/SpacePointFormationError.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <algorithm>
#include <cmath>

namespace Acts {

Result<double> StripSpacePointBuilder::computeClusterPairDistance(
    const Vector3& globalCluster1, const Vector3& globalCluster2,
    const ClusterPairingOptions& options) {
  // Check if measurements are close enough to each other
  if ((globalCluster1 - globalCluster2).norm() > options.maxDistance) {
    return Result<double>::failure(
        SpacePointFormationError::ClusterPairDistanceExceeded);
  }

  const Vector3 vertexToCluster1 = globalCluster1 - options.vertex;
  const Vector3 vertexToCluster2 = globalCluster2 - options.vertex;

  // Calculate the angles of the vectors
  const double phi1 = VectorHelpers::phi(vertexToCluster1);
  const double theta1 = VectorHelpers::theta(vertexToCluster1);
  const double phi2 = VectorHelpers::phi(vertexToCluster2);
  const double theta2 = VectorHelpers::theta(vertexToCluster2);

  // Calculate the squared difference between the theta angles
  const double diffTheta = std::abs(theta1 - theta2);
  if (diffTheta > options.maxAngleTheta) {
    return Result<double>::failure(
        SpacePointFormationError::ClusterPairThetaDistanceExceeded);
  }

  // Calculate the squared difference between the phi angles
  const double diffPhi = std::abs(phi1 - phi2);
  if (diffPhi > options.maxAnglePhi) {
    return Result<double>::failure(
        SpacePointFormationError::ClusterPairPhiDistanceExceeded);
  }

  // Return the squared distance between both vector
  const double distance2 = square(diffTheta) + square(diffPhi);
  return Result<double>::success(distance2);
}

Result<Vector3> StripSpacePointBuilder::computeCosmicSpacePoint(
    const StripEnds& stripEnds1, const StripEnds& stripEnds2,
    const CosmicOptions& options) {
  // This approach assumes that no vertex is available. This option aims to
  // approximate the space points from cosmic data.
  // The underlying assumption is that the best point is given by the closest
  // distance between both lines describing the SDEs. The point x on the first
  // SDE is parametrized as a + lambda0 * q with the top end a of the strip and
  // the vector q = a - (bottom end of the strip). An analogous parametrization
  // is performed of the second SDE with y = c  + lambda1 * r. x get resolved by
  // resolving lambda0 from the condition that |x-y| is the shortest distance
  // between two skew lines.

  // Minimising |x - y|^2 over lambda0 and lambda1 gives
  //   lambda0 * (q.q) - lambda1 * (q.r) = ac.q
  //   lambda0 * (q.r) - lambda1 * (r.r) = ac.r
  // with ac = c - a, which resolves to the lambdas below.

  const Vector3 firstBtmToTop = stripEnds1.top - stripEnds1.bottom;
  const Vector3 secondBtmToTop = stripEnds2.top - stripEnds2.bottom;

  const Vector3 ac = stripEnds2.top - stripEnds1.top;
  const double qq = firstBtmToTop.dot(firstBtmToTop);
  const double rr = secondBtmToTop.dot(secondBtmToTop);
  const double qr = firstBtmToTop.dot(secondBtmToTop);
  const double acq = ac.dot(firstBtmToTop);
  const double acr = ac.dot(secondBtmToTop);

  // By Lagrange's identity this is (q.q) * (r.r) * sin^2(angle between the
  // strips), so the check below is on the opening angle of the two strips and
  // is independent of their length. It is never negative, and zero for a strip
  // of zero length.
  const double denom = qq * rr - qr * qr;
  if (denom <= options.tolerance * qq * rr) {
    return Result<Vector3>::failure(
        SpacePointFormationError::CosmicToleranceNotMet);
  }

  const double lambda0 = (acq * rr - acr * qr) / denom;
  const double lambda1 = (acq * qr - acr * qq) / denom;

  // The crossing has to lie on both strips. This is the stereo matching
  // condition and the only compatibility test available without a vertex.
  // `2 * lambda + 1` centres the on-strip range (-1, 0), so the tolerance is a
  // fraction of the strip length as in the constrained formation.
  const double limit = 1 + options.stripLengthTolerance;
  if (std::abs(2 * lambda0 + 1) > limit || std::abs(2 * lambda1 + 1) > limit) {
    return Result<Vector3>::failure(SpacePointFormationError::OutsideLimits);
  }

  const Vector3 spacePoint = stripEnds1.top + lambda0 * firstBtmToTop;
  return Result<Vector3>::success(spacePoint);
}

StripSpacePointBuilder::ConstrainedStripCache
StripSpacePointBuilder::makeConstrainedStripCache(
    const StripEnds& stripEnds, const ConstrainedOptions& options) {
  ConstrainedStripCache cache;
  cache.btmToTop = stripEnds.top - stripEnds.bottom;
  cache.mid = 0.5 * (stripEnds.top + stripEnds.bottom);
  cache.vtxToMid2 = 2 * (cache.mid - options.vertex);
  cache.normal = cache.btmToTop.cross(cache.vtxToMid2);
  cache.invLength = 1. / cache.btmToTop.norm();
  return cache;
}

std::optional<SpacePointFormationError>
StripSpacePointBuilder::computeConstrainedSpacePoint(
    const ConstrainedStripCache& cache1, const ConstrainedStripCache& cache2,
    const ConstrainedOptions& options, Vector3& spacePoint) {
  /// The following algorithm is meant for finding the position on the first
  /// strip if there is a corresponding Measurement on the second strip. The
  /// resulting point is a point x on the first surfaces. This point is
  /// along a line between the points a (top end of the strip)
  /// and b (bottom end of the strip). The location can be parametrized as
  /// 	2 * x = (1 + m) a + (1 - m) b
  /// as function of the scalar m. m is a parameter in the interval
  /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
  /// from the vertex to the Measurement on the second strip y is needed to be a
  /// multiple k of the vector from vertex to the hit on the first strip x.
  /// As a consequence of this demand y = k * x needs to be on the
  /// connecting line between the top (c) and bottom (d) end of
  /// the second strip. If both measurements correspond to each other, the
  /// condition
  /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
  /// needs to be fulfilled. Inserting the first equation into this
  /// equation leads to the condition for m as given in the following
  /// algorithm and therefore to the calculation of x.
  /// The same calculation can be repeated for y. Its corresponding
  /// parameter will be named n.

  // Regular limit of the absolute values of m and n
  const double limit = 1 + options.stripLengthTolerance;
  // Relaxed limit, allowing the trajectory to be shifted so that a point just
  // beyond the end of a strip is still accepted. The gap tolerance is a length,
  // so it is expressed in units of the strip it is applied to.
  const double relaxedLimit1 =
      limit + options.stripLengthGapTolerance * cache1.invLength;

  // Comparing numerator against denominator rejects the bulk of the candidate
  // pairs without dividing
  const double mNumerator = -cache1.vtxToMid2.dot(cache2.normal);
  const double mDenominator = cache1.btmToTop.dot(cache2.normal);
  if (std::abs(mNumerator) > std::abs(mDenominator) * relaxedLimit1) {
    return SpacePointFormationError::OutsideRelaxedLimits;
  }

  const double relaxedLimit2 =
      limit + options.stripLengthGapTolerance * cache2.invLength;
  const double nNumerator = -cache2.vtxToMid2.dot(cache1.normal);
  const double nDenominator = cache2.btmToTop.dot(cache1.normal);
  if (std::abs(nNumerator) > std::abs(nDenominator) * relaxedLimit2) {
    return SpacePointFormationError::OutsideRelaxedLimits;
  }

  // Parameter that determines the hit position on the first strip
  double m = mNumerator / mDenominator;
  // Parameter that determines the hit position on the second strip
  double n = nNumerator / nDenominator;

  if (std::abs(m) > limit || std::abs(n) > limit) {
    /// The overshoot is the amount by which m or n falls outside (-1, 1). The
    /// worse of the two, compared after projecting n onto the first strip, is
    /// set to +/-1 and the other moved towards 0 by the same shift. Only one
    /// has to be outside: the projection drags the other along.
    /// @note This can be read as a shift of the particle's trajectory, and so
    /// of the vertex, or equivalently as a change of its slope.

    // Calculate the scaling factor to project lengths of the second strip on
    // the first strip
    const double secOnFirstScale = cache1.btmToTop.dot(cache2.btmToTop) *
                                   cache1.invLength * cache1.invLength;

    if (m > limit || n > limit) {
      // Overshoot beyond the top end. Resolve the worse of the two, after
      // projecting the one on the second strip onto the first.
      const double biggerOvershoot = std::max(m - 1, (n - 1) * secOnFirstScale);
      // Move m and n towards 0
      m -= biggerOvershoot;
      n -= biggerOvershoot / secOnFirstScale;
    } else {
      // Overshoot beyond the bottom end
      const double biggerOvershoot =
          std::max(-(m + 1), -(n + 1) * secOnFirstScale);
      m += biggerOvershoot;
      n += biggerOvershoot / secOnFirstScale;
    }

    // Check if this recovered the space point
    if (std::abs(m) > limit || std::abs(n) > limit) {
      return SpacePointFormationError::OutsideLimits;
    }
  }

  spacePoint = cache1.mid + 0.5 * m * cache1.btmToTop;
  return std::nullopt;
}

Result<Vector3> StripSpacePointBuilder::computeConstrainedSpacePoint(
    const ConstrainedStripCache& cache1, const ConstrainedStripCache& cache2,
    const ConstrainedOptions& options) {
  Vector3 spacePoint = Vector3::Zero();
  const std::optional<SpacePointFormationError> error =
      computeConstrainedSpacePoint(cache1, cache2, options, spacePoint);
  if (error.has_value()) {
    return Result<Vector3>::failure(*error);
  }
  return Result<Vector3>::success(spacePoint);
}

Result<Vector3> StripSpacePointBuilder::computeConstrainedSpacePoint(
    const StripEnds& stripEnds1, const StripEnds& stripEnds2,
    const ConstrainedOptions& options) {
  return computeConstrainedSpacePoint(
      makeConstrainedStripCache(stripEnds1, options),
      makeConstrainedStripCache(stripEnds2, options), options);
}

SquareMatrix2 StripSpacePointBuilder::computeCovarianceZR(
    const RotationMatrix3& rotLocalToGlobal, const Vector3& spacePoint,
    const double localCov1, const double localCov2, const double theta) {
  // Invert the information matrix of the two measurements, strip1 along (1, 0)
  // and strip2 along (cos(theta), sin(theta)). Strip2 adds nothing to the
  // precision direction, it only fixes where along strip1 the crossing sits,
  // to 1/sin(theta) of a pitch
  const double sinTheta = std::sin(theta);
  const double cosTheta = std::cos(theta);
  const double varAlongStrip =
      (localCov2 + localCov1 * cosTheta * cosTheta) / (sinTheta * sinTheta);

  // Dropping a correlation of -localCov1 * cot(theta), whose sense an unsigned
  // theta does not carry
  const SquareMatrix2 localCov = Vector2(localCov1, varAlongStrip).asDiagonal();

  return PixelSpacePointBuilder::computeCovarianceZR(rotLocalToGlobal,
                                                     spacePoint, localCov);
}

Vector2 StripSpacePointBuilder::computeVarianceZR(const GeometryContext& gctx,
                                                  const Surface& surface1,
                                                  const Vector3& spacePoint,
                                                  const double localCov1,
                                                  const double localCov2,
                                                  const double theta) {
  return computeCovarianceZR(
             surface1.referenceFrame(gctx, spacePoint, Vector3::Zero()),
             spacePoint, localCov1, localCov2, theta)
      .diagonal();
}

}  // namespace Acts
