// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/SpacePointBuilder.hpp"

namespace Acts {

/// @brief Structure for easier bookkeeping of clusters.
struct SingleHitSpacePoint
{
  /// Storage of the cluster on a surface
  const PlanarModuleCluster* clusterModule;
  /// Storage of a space point.
  Vector3D spacePoint;

  /// @brief Getter of the first element in @p spacePoint
  /// @return First element in @p spacePoint
  double
  x() const
  {
    return spacePoint(0);
  }

  /// @brief Getter of the second element in @p spacePoint
  /// @return Second element in @p spacePoint
  double
  y() const
  {
    return spacePoint(1);
  }

  /// @brief Getter of the third element in @p spacePoint
  /// @return Third element in @p spacePoint
  double
  z() const
  {
    return spacePoint(2);
  }
};

/// @class OneHitSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the clusters pixel detectors need further treatment. This class takes
/// the digitized clusters on a pixel detector element and provides the
/// corresponding space point.
///
template <>
class SpacePointBuilder<SingleHitSpacePoint>
{
public:
  /// Default constructor
  SpacePointBuilder<SingleHitSpacePoint>() = default;

  /// @brief Calculates the space points out of a given collection of clusters
  /// and stores the results
  /// @param cluster vector of clusters
  /// @param spacePointStorage storage of the results
  void
  calculateSpacePoints(
      const std::vector<const PlanarModuleCluster*>& clusters,
      std::vector<SingleHitSpacePoint>&              spacePointStorage) const;

protected:
  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2D
  localCoords(const PlanarModuleCluster& cluster) const;

  /// @brief Getter method for the global coordinates of a cluster
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the global coordinates of the cluster
  Vector3D
  globalCoords(const PlanarModuleCluster& cluster) const;
};

}  // namespace Acts
