// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/Utilities/Result.hpp"

#include <span>

namespace Acts {

class GeometryContext;
class Surface;

namespace StripSpacePointBuilder {

/// @brief Collection of cluster pairing options
struct ClusterPairingOptions final {
  /// vertex position
  Vector3 vertex = Vector3::Zero();
  /// Accepted distance between two clusters
  double maxDistance = 100 * UnitConstants::mm;
  /// Accepted absolute difference in theta for two clusters
  double maxAngleTheta = 1 * UnitConstants::rad;
  /// Accepted absolute difference in phi for two clusters
  double maxAnglePhi = 1 * UnitConstants::rad;
};

/// @brief Collection of cosmic space point options
struct CosmicOptions final {
  /// Numerical tolerance for the calculation
  double tolerance = 1e-6;
};

/// @brief Collection of constrained space point options
struct ConstrainedOptions final {
  /// Position of the vertex
  Vector3 vertex = Vector3::Zero();
  /// Tolerance scaling factor on the strip detector element length
  double stripLengthTolerance = 0.01;
  /// Tolerance scaling factor of the gap between strip detector elements
  double stripLengthGapTolerance = 0.01;
};

/// @brief Strip cluster details
struct StripEnds final {
  /// Top end of the strip cluster
  Vector3 top = Vector3::Zero();
  /// Bottom end of the strip cluster
  Vector3 bottom = Vector3::Zero();
};

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two measurements
///
/// @param globalCluster1 Global position of the measurements on the first strip
/// @param globalCluster2 Global position of the measurements on the second strip
/// @param options Pairing options
///
/// @return If available, squared sum within configuration parameters
Result<double> computeClusterPairDistance(const Vector3& globalCluster1,
                                          const Vector3& globalCluster2,
                                          const ClusterPairingOptions& options);

/// @param stripEnds1 The ends of first strip
/// @param stripEnds2 The ends of second strip
/// @param options The cosmic options
///
/// @return If available, the calculated space point
Result<Vector3> computeCosmicSpacePoint(const StripEnds& stripEnds1,
                                        const StripEnds& stripEnds2,
                                        const CosmicOptions& options);

/// @param stripEnds1 The ends of first strip
/// @param stripEnds2 The ends of second strip
/// @param options The constrained options
///
/// @return If available, the calculated space point
Result<Vector3> computeConstrainedSpacePoint(const StripEnds& stripEnds1,
                                             const StripEnds& stripEnds2,
                                             const ConstrainedOptions& options);

/// @brief Calculate the z and r covariance from the front and back SourceLink in the strip SP formation
///
/// @param gctx The current geometry context object, e.g. alignment
/// @param surface1 The surface of the first strip
/// @param spacePoint The space point
/// @param localCov1 Local covariance of the first strip
/// @param localCov2 Local covariance of the second strip
/// @param theta The angle between the two strips
///
/// @return (z, r) components of the global covariance
Vector2 computeVarianceZR(const GeometryContext& gctx, const Surface& surface1,
                          const Vector3& spacePoint, double localCov1,
                          double localCov2, double theta);

inline StripSpacePointCalibrationDetailsDerived
deriveStripSpacePointCalibrationDetails(
    const StripSpacePointCalibrationDetails& sp) {
  const auto& ohv = sp.outerStripHalfVector;
  const auto& ihv = sp.innerStripHalfVector;
  const auto& scd = sp.stripSeparation;

  return StripSpacePointCalibrationDetailsDerived{
      .stripSeparationCrossInnerHalfVector = {scd[1] * ihv[2] - scd[2] * ihv[1],
                                              scd[2] * ihv[0] - scd[0] * ihv[2],
                                              scd[0] * ohv[1] -
                                                  scd[1] * ohv[0]},
      .stripSeparationCrossOuterHalfVector = {scd[1] * ohv[2] - scd[2] * ohv[1],
                                              scd[2] * ohv[0] - scd[0] * ohv[2],
                                              scd[0] * ohv[1] -
                                                  scd[1] * ohv[0]},
      .innerCrossOuterStripHalfVector = {ihv[1] * ohv[2] - ihv[2] * ohv[1],
                                         ihv[2] * ohv[0] - ihv[0] * ohv[2],
                                         ihv[0] * ohv[1] - ihv[1] * ohv[0]},
      .outerStripCenter = sp.outerStripCenter,
      .outerStripHalfVector = ohv,
  };
}

inline bool calibrateStripSpacePoint(
    const StripSpacePointCalibrationDetails& sp,
    std::span<const float, 3> direction, std::span<float, 3> calibrated,
    float tolerance) {
  const auto& scd = sp.stripSeparation;
  const auto& ohv = sp.outerStripHalfVector;
  const auto& ihv = sp.innerStripHalfVector;

  // dOuter = outerStripHalfVector cross direction (reused for both scale and
  // sInner)
  const std::array<float, 3> dOuter = {
      ohv[1] * direction[2] - ohv[2] * direction[1],
      ohv[2] * direction[0] - ohv[0] * direction[2],
      ohv[0] * direction[1] - ohv[1] * direction[0]};

  // scale = innerStripHalfVector dot d1
  const float scale =
      ihv[0] * dOuter[0] + ihv[1] * dOuter[1] + ihv[2] * dOuter[2];

  // sInner = stripSeparation dot dOuter
  // Check if sInner is inside the inner detector element
  // TODO should this be `sOuter`?
  const float sInner =
      scd[0] * dOuter[0] + scd[1] * dOuter[1] + scd[2] * dOuter[2];
  if (std::abs(sInner) > std::abs(scale) * tolerance) {
    return false;
  }

  // dInner = innerStripHalfVector cross direction (only computed if check 1
  // passed)
  const std::array<float, 3> dInner = {
      ihv[1] * direction[2] - ihv[2] * direction[1],
      ihv[2] * direction[0] - ihv[0] * direction[2],
      ihv[0] * direction[1] - ihv[1] * direction[0]};

  // sOuter = stripSeparation dot dInner
  // Check if sOuter is inside the outer detector element
  // TODO should this be `sInner`?
  const float sOuter =
      scd[0] * dInner[0] + scd[1] * dInner[1] + scd[2] * dInner[2];
  if (std::abs(sOuter) > std::abs(scale) * tolerance) {
    return false;
  }

  // Corrected position using the outer strip center and direction
  // TODO use inner?
  const auto& oc = sp.outerStripCenter;
  const float sOuterNorm = sOuter / scale;
  calibrated[0] = oc[0] + ohv[0] * sOuterNorm;
  calibrated[1] = oc[1] + ohv[1] * sOuterNorm;
  calibrated[2] = oc[2] + ohv[2] * sOuterNorm;
  return true;
}

inline bool calibrateStripSpacePoint(
    const StripSpacePointCalibrationDetailsDerived& sp,
    std::span<const float, 3> direction, std::span<float, 3> calibrated,
    float tolerance) {
  const auto& oc = sp.outerStripCenter;
  const auto& ohv = sp.outerStripHalfVector;
  const auto& sCrossIhv = sp.stripSeparationCrossInnerHalfVector;
  const auto& sCrossOhv = sp.stripSeparationCrossOuterHalfVector;
  const auto& ihvCrossOhv = sp.innerCrossOuterStripHalfVector;

  // scale = innerStripHalfVector dot (outerStripHalfVector cross direction)
  const float scale = direction[0] * ihvCrossOhv[0] +
                      direction[1] * ihvCrossOhv[1] +
                      direction[2] * ihvCrossOhv[2];

  // sInner = stripSeparation dot (outerStripHalfVector cross direction)
  // Check if direction is inside the inner detector element
  // TODO should this rather use `sCrossIhv`?
  const float sInner = direction[0] * sCrossOhv[0] +
                       direction[1] * sCrossOhv[1] +
                       direction[2] * sCrossOhv[2];
  if (std::abs(sInner) > std::abs(scale) * tolerance) {
    return false;
  }

  // sOuter = stripSeparation dot (innerStripHalfVector cross direction)
  // Check if direction is inside the outer detector element
  // TODO should this rather use `sCrossOhv`?
  const float sOuter = direction[0] * sCrossIhv[0] +
                       direction[1] * sCrossIhv[1] +
                       direction[2] * sCrossIhv[2];
  if (std::abs(sOuter) > std::abs(scale) * tolerance) {
    return false;
  }

  // Corrected position using the outer strip center and direction
  // TODO use inner?
  const float sOuterNorm = sOuter / scale;
  calibrated[0] = oc[0] + ohv[0] * sOuterNorm;
  calibrated[1] = oc[1] + ohv[1] * sOuterNorm;
  calibrated[2] = oc[2] + ohv[2] * sOuterNorm;
  return true;
}

}  // namespace StripSpacePointBuilder

}  // namespace Acts
