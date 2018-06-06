// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"
#include "Acts/Tools/ISpacePointBuilder.hpp"

namespace Acts {

/// @brief Structure for easier bookkeeping of clusters.
struct SingleHitSpacePoint
{
  /// Storage of the cluster on a surface
  PlanarModuleCluster const* clusterModule;
  /// Storage of a space point. Zero vector indicates unset point
  Vector3D spacePoint = {0., 0., 0.};

  /// @brief Getter of the first element in @p spacePoint
  /// @return First element in @p spacePoint
  double
  x()
  {
    return spacePoint(0);
  }

  /// @brief Getter of the second element in @p spacePoint
  /// @return Second element in @p spacePoint
  double
  y()
  {
    return spacePoint(1);
  }

  /// @brief Getter of the third element in @p spacePoint
  /// @return Third element in @p spacePoint
  double
  z()
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
  SpacePointBuilder<SingleHitSpacePoint>();

  /// @brief Adds clusters on surfaces and stores them
  /// @param spacePointStorage storage of clusters and the therewith resulting
  /// space points
  /// @param clusters vector of clusters on the surface
  void
  addClusters(std::vector<SingleHitSpacePoint>& spacePointStorage,
              const std::vector<Acts::PlanarModuleCluster const*>& clusters);

  /// @brief Calculates the space points out of a given collection of clusters
  /// and stores the results
  /// @param spacePointStorage storage of the data
  void
  calculateSpacePoints(std::vector<SingleHitSpacePoint>& spacePointStorage);

protected:
  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2D
  localCoords(const PlanarModuleCluster& cluster);

  /// @brief Getter method for the global coordinates of a cluster
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the global coordinates of the cluster
  Vector3D
  globalCoords(const PlanarModuleCluster& cluster);
};

}  // namespace Acts
