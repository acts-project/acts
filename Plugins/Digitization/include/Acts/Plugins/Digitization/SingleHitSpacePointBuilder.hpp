// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/SpacePointBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {

/// @class OneHitSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the clusters pixel detectors need further treatment. This class takes
/// the digitized clusters on a pixel detector element and provides the
/// corresponding space point.
///
template <typename Cluster>
class SpacePointBuilder<SpacePoint<Cluster>> {
 public:
  /// Default constructor
  SpacePointBuilder<SpacePoint<Cluster>>() = default;

  /// @brief Calculates the space points out of a given collection of clusters
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param cluster vector of clusters
  /// @param spacePointStorage storage of the results
  void calculateSpacePoints(
      const GeometryContext& gctx, const std::vector<const Cluster*>& clusters,
      std::vector<SpacePoint<Cluster>>& spacePointStorage) const;

 protected:
  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  ///
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2D localCoords(const Cluster& cluster) const;

  /// @brief Getter method for the global coordinates of a cluster
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the global coordinates of the cluster
  Vector3D globalCoords(const GeometryContext& gctx,
                                const Cluster& cluster) const;
};
}  // namespace Acts
#include "Acts/Plugins/Digitization/detail/SingleHitSpacePointBuilder.ipp"