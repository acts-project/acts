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

/// @brief Collection of cosmic spacepoint options
struct CosmicOptions final {
  /// Numerical tolerance for the calculation
  double tolerance = 1e-6;
};

/// @brief Collection of constrained spacepoint options
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
/// @return If available, the calculated spacepoint
Result<Vector3> computeCosmicSpacePoint(const StripEnds& stripEnds1,
                                        const StripEnds& stripEnds2,
                                        const CosmicOptions& options);

/// @param stripEnds1 The ends of first strip
/// @param stripEnds2 The ends of second strip
/// @param options The constrained options
///
/// @return If available, the calculated spacepoint
Result<Vector3> computeConstrainedSpacePoint(const StripEnds& stripEnds1,
                                             const StripEnds& stripEnds2,
                                             const ConstrainedOptions& options);

/// @brief Calculate the z and r covariance from the front and back SourceLink in the strip SP formation
///
/// @param gctx The current geometry context object, e.g. alignment
/// @param surface1 The surface of the first strip
/// @param spacePoint The spacepoint
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
