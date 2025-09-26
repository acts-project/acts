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

namespace Acts {

class GeometryContext;
class Surface;

namespace StripSpacePointBuilder {

struct ClusterPairingOptions {
  /// vertex position
  Vector3 vertex = Vector3::Zero();
  /// Accepted distance between two clusters
  double maxDistance = 100 * UnitConstants::mm;
  /// Accepted squared difference in theta for two clusters
  double maxAngleTheta2 = 1;
  /// Accepted squared difference in phi for two clusters
  double maxAnglePhi2 = 1;
};

struct StripEnds {
  Vector3 top = Vector3::Zero();
  Vector3 bottom = Vector3::Zero();
};

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two measurements
///
/// @param globalCluster1 Global position of the first measurements
/// @param globalCluster2 Global position the second measurements
/// @param options Pairing options
///
/// @return Result with the squared sum within configuration parameters.
std::optional<double> computeClusterPairDistance(
    const Vector3& globalCluster1, const Vector3& globalCluster2,
    const ClusterPairingOptions& options);

/// @param stripEnds1 The ends of one strip
/// @param stripEnds2 The ends of another strip
/// @param spacePoint The calculated space point
/// @return whether the space point calculation was successful
bool computeCosmicSpacePoint(const StripEnds& stripEnds1,
                             const StripEnds& stripEnds2, Vector3& spacePoint);

/// @param stripEnds1 The ends of one strip
/// @param stripEnds2 The ends of another strip
/// @param vertex Position of the vertex
/// @param spacePoint The calculated space point
/// @param stripLengthTolerance Tolerance scaling factor on the
/// strip detector element length
/// @param stripLengthGapTolerance Tolerance scaling factor of
/// the gap between strip detector elements
/// @return whether the space point calculation was successful
bool computeConstrainedSpacePoint(const StripEnds& stripEnds1,
                                  const StripEnds& stripEnds2,
                                  const Vector3& vertex, Vector3& spacePoint,
                                  double stripLengthTolerance,
                                  double stripLengthGapTolerance);

/// @brief Calculate the z and r covariance from the front and back SourceLink in the strip SP formation
/// @param gctx The current geometry context object, e.g. alignment
/// @param surface1 The surface of the first strip
/// @param spacePoint The space point
/// @param localCov1 Local covariance of the first strip
/// @param localCov2 Local covariance of the second strip
/// @param theta The angle between the two strips
/// @return (z, r) components of the global covariance
Vector2 computeVarianceZR(const GeometryContext& gctx, const Surface& surface1,
                          const Vector3& spacePoint, double localCov1,
                          double localCov2, double theta);

}  // namespace StripSpacePointBuilder

}  // namespace Acts
