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
#include "Acts/SpacePointFormation/SpacePointFormationError.hpp"
#include "Acts/Utilities/Result.hpp"

#include <optional>

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
  /// Minimum squared sine of the angle between the two strips. Pairs below this
  /// are too close to parallel for the closest approach to be well defined.
  double tolerance = 1e-6;
  /// Tolerance scaling factor on the strip detector element length
  double stripLengthTolerance = 0.01;
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

/// @brief Calculates the space point from the closest approach of two strips
///
/// Rejects the pair if the closest approach does not lie on both strips.
///
/// @param stripEnds1 The ends of first strip
/// @param stripEnds2 The ends of second strip
/// @param options The cosmic options
///
/// @return If available, the calculated space point
Result<Vector3> computeCosmicSpacePoint(const StripEnds& stripEnds1,
                                        const StripEnds& stripEnds2,
                                        const CosmicOptions& options);

/// @brief Per-strip cache for the constrained space point formation
///
/// Depends on one strip only, so it can be filled once per strip and reused
/// for every candidate partner on the opposite side.
struct ConstrainedStripCache final {
  /// Vector pointing from the bottom to the top end of the strip
  Vector3 btmToTop = Vector3::Zero();
  /// Twice the vector pointing from the vertex to the midpoint of the strip
  Vector3 vtxToMid2 = Vector3::Zero();
  /// Normal of the plane spanned by the strip and the trajectory,
  /// `btmToTop x vtxToMid2`
  Vector3 normal = Vector3::Zero();
  /// Midpoint of the strip
  Vector3 mid = Vector3::Zero();
  /// Inverse of the strip length
  double invLength = 0;
};

/// @brief Precompute the per-strip quantities of the constrained formation
///
/// @param stripEnds The ends of the strip
/// @param options The constrained options; only the vertex is used
///
/// @return The cache to hand to computeConstrainedSpacePoint
ConstrainedStripCache makeConstrainedStripCache(
    const StripEnds& stripEnds, const ConstrainedOptions& options);

/// @param stripEnds1 The ends of first strip
/// @param stripEnds2 The ends of second strip
/// @param options The constrained options
///
/// @return If available, the calculated space point
Result<Vector3> computeConstrainedSpacePoint(const StripEnds& stripEnds1,
                                             const StripEnds& stripEnds2,
                                             const ConstrainedOptions& options);

/// @brief Space point of a strip pair, from precomputed per-strip caches
///
/// Lets the caller hoist the per-strip work out of the loop over candidate
/// pairs. Both caches must have been filled with the same options;
/// `options.vertex` is not read again here, it is already encoded in them.
///
/// @param cache1 The cache of the first strip
/// @param cache2 The cache of the second strip
/// @param options The constrained options
/// @param spacePoint Set to the calculated space point, untouched on failure
///
/// @return Why the pair was rejected, or nullopt if a space point was formed
std::optional<SpacePointFormationError> computeConstrainedSpacePoint(
    const ConstrainedStripCache& cache1, const ConstrainedStripCache& cache2,
    const ConstrainedOptions& options, Vector3& spacePoint);

/// @brief Space point of a strip pair, from precomputed per-strip caches
///
/// As above, packed into a Result.
///
/// @param cache1 The cache of the first strip
/// @param cache2 The cache of the second strip
/// @param options The constrained options
///
/// @return If available, the calculated space point
Result<Vector3> computeConstrainedSpacePoint(
    const ConstrainedStripCache& cache1, const ConstrainedStripCache& cache2,
    const ConstrainedOptions& options);

/// @brief Transform the crossing covariance of a strip pair to z and r
///
/// See Acts::PixelSpacePointBuilder::computeCovarianceZR for the requirement
/// on the reference frame
///
/// @param rotLocalToGlobal The reference frame of the first strip's surface
/// @param spacePoint The space point
/// @param localCov1 Local covariance of the first strip
/// @param localCov2 Local covariance of the second strip
/// @param theta The angle between the two strips
///
/// @return The (z, r) covariance
SquareMatrix2 computeCovarianceZR(const RotationMatrix3& rotLocalToGlobal,
                                  const Vector3& spacePoint, double localCov1,
                                  double localCov2, double theta);

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
/// @deprecated Use computeCovarianceZR with a reference frame from
///             Surface::referenceFrame
[[deprecated(
    "Use computeCovarianceZR with a reference frame from "
    "Surface::referenceFrame")]]
Vector2 computeVarianceZR(const GeometryContext& gctx, const Surface& surface1,
                          const Vector3& spacePoint, double localCov1,
                          double localCov2, double theta);

}  // namespace StripSpacePointBuilder

}  // namespace Acts
