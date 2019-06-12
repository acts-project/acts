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
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// @class TwoHitsSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the hits strip detectors need further treatment. This class takes
/// the digitized clusters and combines them on two different detector elements
/// to a result of the combined detector element. The class is intended to
/// handle strip detector elements in particular.
///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///
template <typename Cluster>
class SpacePointBuilder<SpacePoint<Cluster>> {
 public:
  /// @brief Configuration of the class to steer its behaviour
  struct DoubleHitSpacePointConfig {
    /// Accepted difference in eta for two clusters
    double diffTheta2 = 1.;
    /// Accepted difference in phi for two clusters
    double diffPhi2 = 1.;
    /// Accepted distance between two clusters
    double diffDist = 100. * UnitConstants::mm;
    /// Allowed increase of strip length
    double stripLengthTolerance = 0.01;
    /// Allowed increase of strip length wrt gaps between strips
    double stripLengthGapTolerance = 0.01;
    /// Assumed position of the vertex
    Vector3D vertex = {0., 0., 0.};
    /// Perform the perpendicular projection for space point finding
    bool usePerpProj = false;
  };

  /// Constructor
  /// @param cfg Specific config that will be used instead of the default values
  SpacePointBuilder(DoubleHitSpacePointConfig cfg);

  /// @brief Searches possible combinations of two clusters on different
  /// surfaces that may come from the same particles
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clustersFront vector of clusters on a surface
  /// @param clustersBack vector of clusters on another surface
  /// @param clusterPairs storage of the cluster pairs
  /// @note The structure of @p clustersFront and @p clustersBack is meant to be
  /// clusters[Independent clusters on a single surface]
  void makeClusterPairs(const GeometryContext& gctx,
                        const std::vector<const Cluster*>& clustersFront,
                        const std::vector<const Cluster*>& clustersBack,
                        std::vector<std::pair<const Cluster*, const Cluster*>>&
                            clusterPairs) const;

  /// @brief Calculates the space points out of a given collection of clusters
  /// on several strip detectors and stores the data
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clusterPairs pairs of clusters that are space point candidates
  /// @param spacePoints storage of the results
  /// @note If no configuration is set, the default values will be used
  void calculateSpacePoints(
      const GeometryContext& gctx,
      const std::vector<std::pair<const Cluster*, const Cluster*>>&
          clusterPairs,
      std::vector<SpacePoint<Cluster>>& spacePoints) const;

 private:
  /// Config
  DoubleHitSpacePointConfig m_cfg;

  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2D localCoords(const Cluster& cluster) const;

  /// @brief Getter method for the global coordinates of a cluster
  /// @param cluster object related to the cluster that holds the necessary
  /// information
  /// @return vector of the global coordinates of the cluster
  Vector3D globalCoords(const GeometryContext& gctx,
                        const Cluster& cluster) const;

  /// @brief Calculates the top and bottom ends of a SDE
  /// that corresponds to a given hit
  /// @param cluster object that stores the information about the hit
  /// @return vectors to the top and bottom end of the SDE
  std::pair<Vector3D, Vector3D> endsOfStrip(const GeometryContext& gctx,
                                            const Cluster& cluster) const;
};
#include "Acts/Plugins/Digitization/detail/DoubleHitSpacePointBuilder.ipp"
}  // namespace Acts