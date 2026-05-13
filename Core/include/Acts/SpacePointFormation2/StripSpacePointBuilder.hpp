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
#include "Acts/Utilities/Result.hpp"

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
///
/// The constrained space-point algorithm parametrises the hit position on
/// each strip by a scalar:
///   strip 1: 2 * x_1 = (1 + m) * top_1 + (1 - m) * bottom_1
///   strip 2: 2 * x_2 = (1 + n) * top_2 + (1 - n) * bottom_2
/// `m = n = 0` is the strip midpoint; `|m| = 1` (resp. `|n| = 1`) is the
/// strip end. A valid hit requires `|m|, |n| <= limit`.
struct ConstrainedOptions final {
  /// Position of the vertex
  Vector3 vertex = Vector3::Zero();
  /// Additive slack on the `|m|, |n| <= 1` acceptance: limit = 1 + tol.
  double stripLengthTolerance = 0.01;
  /// Additional slack used inside the recovery branch when `|m|` or `|n|`
  /// is just past the limit.
  double stripLengthGapTolerance = 0.01;
  /// Additional acceptance allowance that scales with the geometric gap
  /// between the two stereo wafer faces. When > 0, the limit becomes
  ///   1 + stripLengthTolerance
  ///     + stripGapParameter * |gapVec| / sin(stereo) / stripHalfLength,
  /// where `gapVec` is the vector connecting the two strip midpoints,
  /// `stereo` is the angle between the two strip directions, and
  /// `stripHalfLength` is half of `|top_1 - bottom_1|`. This is the
  /// geometric-allowance piece of the ATLAS Athena strip pair algorithm
  /// (see `SiSpacePointMakerTool::makeSCT_SpacePoint`'s `m_SCTgapParameter`
  /// usage); it lets high-`|eta|` tracks that physically cross the
  /// wafer-thickness gap survive the acceptance gate. 0 disables it
  /// (the strict `|m|, |n| <= 1 + stripLengthTolerance` cut applies).
  double stripGapParameter = 0.0;
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

}  // namespace StripSpacePointBuilder

}  // namespace Acts
