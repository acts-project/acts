// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation2/StripSpacePointBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/SpacePointFormation2/PixelSpacePointBuilder.hpp"
#include "Acts/SpacePointFormation2/SpacePointFormationError.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

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
  // approximate the spacepoints from cosmic data.
  // The underlying assumption is that the best point is given by the closest
  // distance between both lines describing the SDEs. The point x on the first
  // SDE is parametrized as a + lambda0 * q with the top end a of the strip and
  // the vector q = a - (bottom end of the strip). An analogous parametrization
  // is performed of the second SDE with y = c  + lambda1 * r. x get resolved by
  // resolving lambda0 from the condition that |x-y| is the shortest distance
  // between two skew lines.

  const Vector3 firstBtmToTop = stripEnds1.top - stripEnds1.bottom;
  const Vector3 secondBtmToTop = stripEnds2.top - stripEnds2.bottom;

  const Vector3 ac = stripEnds2.top - stripEnds1.top;
  const double qr = firstBtmToTop.dot(secondBtmToTop);
  const double denom = firstBtmToTop.dot(firstBtmToTop) - qr * qr;
  // Check for numerical stability
  if (std::abs(denom) < options.tolerance) {
    return Result<Vector3>::failure(
        SpacePointFormationError::CosmicToleranceNotMet);
  }

  const double lambda0 =
      (ac.dot(secondBtmToTop) * qr -
       ac.dot(firstBtmToTop) * secondBtmToTop.dot(secondBtmToTop)) /
      denom;
  const Vector3 spacePoint = stripEnds1.top + lambda0 * firstBtmToTop;
  return Result<Vector3>::success(spacePoint);
}

namespace {

/// @brief Struct for variables related to the calculation of spacepoints
struct FormationState {
  /// Vector pointing from bottom to top end of first SDE
  Vector3 firstBtmToTop = Vector3::Zero();
  /// Vector pointing from bottom to top end of second SDE
  Vector3 secondBtmToTop = Vector3::Zero();
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3 vtxToFirstMid2 = Vector3::Zero();
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3 vtxToSecondMid2 = Vector3::Zero();
  /// Cross product between firstBtmToTop and vtxToFirstMid2
  Vector3 firstBtmToTopXvtxToFirstMid2 = Vector3::Zero();
  /// Cross product between secondBtmToTop and vtxToSecondMid2
  Vector3 secondBtmToTopXvtxToSecondMid2 = Vector3::Zero();
  /// Parameter that determines the hit position on the first SDE
  double m = 0;
  /// Parameter that determines the hit position on the second SDE
  double n = 0;
  /// Regular limit of the absolute values of `m` and `n`
  double limit = 1;
};

/// @brief This function performs a straight forward calculation of a space
/// point and returns whether it was successful or not.
///
/// @param stripEnds1 Top and bottom end of the first strip
/// @param stripEnds2 Top and bottom end of the second strip
/// @param vertex Position of the vertex
/// @param spParams Data container of the calculations
/// @param stripLengthTolerance Tolerance scaling factor on the strip detector element length
///
/// @return whether the spacepoint calculation was successful
Result<void> computeConstrainedFormationState(
    const StripSpacePointBuilder::StripEnds& stripEnds1,
    const StripSpacePointBuilder::StripEnds& stripEnds2, const Vector3& vertex,
    FormationState& state, const double stripLengthTolerance) {
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

  state.firstBtmToTop = stripEnds1.top - stripEnds1.bottom;
  state.secondBtmToTop = stripEnds2.top - stripEnds2.bottom;
  state.vtxToFirstMid2 = stripEnds1.top + stripEnds1.bottom - 2 * vertex;
  state.vtxToSecondMid2 = stripEnds2.top + stripEnds2.bottom - 2 * vertex;
  state.firstBtmToTopXvtxToFirstMid2 =
      state.firstBtmToTop.cross(state.vtxToFirstMid2);
  state.secondBtmToTopXvtxToSecondMid2 =
      state.secondBtmToTop.cross(state.vtxToSecondMid2);
  state.m = -state.vtxToFirstMid2.dot(state.secondBtmToTopXvtxToSecondMid2) /
            state.firstBtmToTop.dot(state.secondBtmToTopXvtxToSecondMid2);
  state.n = -state.vtxToSecondMid2.dot(state.firstBtmToTopXvtxToFirstMid2) /
            state.secondBtmToTop.dot(state.firstBtmToTopXvtxToFirstMid2);

  // Set the limit for the parameter
  state.limit = 1 + stripLengthTolerance;

  // Check if m and n can be resolved in the interval (-1, 1)
  if (std::abs(state.m) > state.limit || std::abs(state.n) > state.limit) {
    return Result<void>::failure(SpacePointFormationError::OutsideLimits);
  }
  return Result<void>::success();
}

/// @brief This function tests if a spacepoint can be estimated by a more
/// tolerant treatment of construction. In fact, this function indirectly
/// allows shifts of the vertex.
///
/// @param state container that stores geometric parameters and rules of the spacepoint formation
/// @param stripLengthGapTolerance Tolerance scaling factor of the gap between strip detector elements
///
/// @return whether the spacepoint calculation was successful
Result<void> recoverConstrainedFormationState(
    FormationState& state, const double stripLengthGapTolerance) {
  const double magFirstBtmToTop = state.firstBtmToTop.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  const double relaxedLimit =
      state.limit + stripLengthGapTolerance / magFirstBtmToTop;

  // Check if m is just slightly outside
  if (std::abs(state.m) > relaxedLimit) {
    return Result<void>::failure(
        SpacePointFormationError::OutsideRelaxedLimits);
  }
  // Calculate n if not performed previously
  if (state.n == 0) {
    state.n = -state.vtxToSecondMid2.dot(state.firstBtmToTopXvtxToFirstMid2) /
              state.secondBtmToTop.dot(state.firstBtmToTopXvtxToFirstMid2);
  }
  // Check if n is just slightly outside
  if (std::abs(state.n) > relaxedLimit) {
    return Result<void>::failure(
        SpacePointFormationError::OutsideRelaxedLimits);
  }

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the spacepoint can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other measurement, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory.
  /// This would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  const double secOnFirstScale =
      state.firstBtmToTop.dot(state.secondBtmToTop) / square(magFirstBtmToTop);

  // Check if both overshoots are in the same direction
  if (state.m > 1 && state.n > 1) {
    // Calculate the overshoots
    const double mOvershoot = state.m - 1;
    // Perform projection
    const double nOvershoot = (state.n - 1) * secOnFirstScale;
    // Resolve worse overshoot
    const double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    state.m -= biggerOvershoot;
    state.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the spacepoint
    if (std::abs(state.m) > state.limit || std::abs(state.n) > state.limit) {
      return Result<void>::failure(SpacePointFormationError::OutsideLimits);
    }
    return Result<void>::success();
  }

  // Check if both overshoots are in the same direction
  if (state.m < -1 && state.n < -1) {
    // Calculate the overshoots
    const double mOvershoot = -(state.m + 1);
    // Perform projection
    const double nOvershoot = -(state.n + 1) * secOnFirstScale;
    // Resolve worse overshoot
    const double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    state.m += biggerOvershoot;
    state.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the spacepoint
    if (std::abs(state.m) > state.limit || std::abs(state.n) > state.limit) {
      return Result<void>::failure(SpacePointFormationError::OutsideLimits);
    }
    return Result<void>::success();
  }

  // No solution could be found
  return Result<void>::failure(SpacePointFormationError::NoSolutionFound);
}

}  // namespace

Result<Vector3> StripSpacePointBuilder::computeConstrainedSpacePoint(
    const StripEnds& stripEnds1, const StripEnds& stripEnds2,
    const ConstrainedOptions& options) {
  FormationState state;

  Result<void> result =
      computeConstrainedFormationState(stripEnds1, stripEnds2, options.vertex,
                                       state, options.stripLengthTolerance);

  if (!result.ok()) {
    result = recoverConstrainedFormationState(state,
                                              options.stripLengthGapTolerance);
  }

  if (!result.ok()) {
    return Result<Vector3>::failure(result.error());
  }

  const Vector3 spacePoint = 0.5 * (stripEnds1.top + stripEnds1.bottom +
                                    state.m * state.firstBtmToTop);
  return Result<Vector3>::success(spacePoint);
}

Vector2 StripSpacePointBuilder::computeVarianceZR(const GeometryContext& gctx,
                                                  const Surface& surface1,
                                                  const Vector3& spacePoint,
                                                  const double localCov1,
                                                  const double localCov2,
                                                  const double theta) {
  const double sinThetaHalf = std::sin(0.5 * theta);
  const double cosThetaHalf = std::cos(0.5 * theta);

  // strip1 and strip2 are tilted at +/- theta/2
  const double var = fastHypot(localCov1, localCov2);
  const double varX = var / (2 * sinThetaHalf);
  const double varY = var / (2 * cosThetaHalf);

  // projection to the surface with strip1.
  const double varX1 = varX * cosThetaHalf + varY * sinThetaHalf;
  const double varY1 = varY * cosThetaHalf + varX * sinThetaHalf;
  const SquareMatrix2 localCov = Vector2(varX1, varY1).asDiagonal();

  return PixelSpacePointBuilder::computeVarianceZR(gctx, surface1, spacePoint,
                                                   localCov);
}

}  // namespace Acts
