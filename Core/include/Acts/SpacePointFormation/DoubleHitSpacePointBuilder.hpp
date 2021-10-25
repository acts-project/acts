// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

namespace Acts {

/// @brief Configuration of the class to steer its behaviour

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

template <typename spacepoint_t, typename cluster_t>
class DoubleHitSpacePointBuilder {
 public:
  /// Constructor
  /// @param cfg Specific config that will be used instead of the default values
  DoubleHitSpacePointBuilder(DoubleHitSpacePointBuilderConfig cfg);
  ///// Default constructor
  DoubleHitSpacePointBuilder() = default;
  /// @brief Searches possible combinations of two clusters on different
  /// surfaces that may come from the same particles
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clustersFront vector of clusters on a surface
  /// @param clustersBack vector of clusters on another surface
  /// @param clusterPairs storage of the cluster pairs
  /// @note The structure of @p clustersFront and @p clustersBack is
  /// meant to be clusters[Independent clusters on a single surface]
  void makeClusterPairs(
      const GeometryContext& gctx,
      const std::vector<const cluster_t*>& clustersFront,
      const std::vector<const cluster_t*>& clustersBack,
      std::vector<std::pair<const cluster_t*, const cluster_t*>>& clusterPairs)
      const;

  /// @brief Calculates the space points out of a given collection of
  /// clusters on several strip detectors and stores the data
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clusterPairs pairs of clusters that are space point
  /// candidates
  /// @param spacePoints storage of the results
  void calculateSpacePoints(
      const GeometryContext& gctx,
      const std::vector<std::pair<const cluster_t*, const cluster_t*>>&
          clusterPairs,
      std::vector<spacepoint_t>& spacePoints) const;

 private:
  /// Config
  DoubleHitSpacePointBuilderConfig m_cfg;

  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  /// @param cluster object related to the cluster that holds the necessary
  ///                information
  /// @return vector of the local coordinates and covariance matrix of the cluster on the surface
  std::pair<Vector2, SymMatrix2> localCoords(const cluster_t& cluster) const;

  /// @brief Getter method for the global coordinates of a cluster
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param cluster object related to the cluster that holds the
  /// necessary information
  /// @return vector of the global coordinates of the cluster
  Vector3 globalPos(const GeometryContext& gctx,
                    const cluster_t& cluster) const;

  /// @brief Calculates the top and bottom ends of a SDE
  /// that corresponds to a given hit
  /// @param gctx The geometry context to use
  /// @param cluster object that stores the information about the hit
  /// @return vectors to the top and bottom end of the SDE
  std::pair<Vector3, Vector3> endsOfStrip(const GeometryContext& gctx,
                                          const cluster_t& cluster) const;

  /// @brief Calculates the global covariance in R and Z directions
  /// @param gctx The geometry context to use
  /// @param geoId The geometry ID of the surface
  /// @param localPos Local coordinates
  /// @param localCov local covariant matrix that is calculated in this function
  /// @return global variances in R and Z directions.
  Vector2 globalCov(const GeometryContext& gctx,
                    const GeometryIdentifier& geoId, const Vector2& localPos,
                    const SymMatrix2& localCov) const;

  /// @brief Gets measurement ID of the cluster
  /// @param cluster object that stores the information about the hit
  /// @return measurement ID
  size_t getMeasurementId(const cluster_t& cluster) const;

  /// @brief Gets local variance of the strip cluster
  /// @param cluster object that sotres the information about the strip hit
  /// @return lcoal variance of the strip
  double getLocVar(const cluster_t& cluster) const;

  /// @brief gets global variations in R and Z directions from two strip clusters.
  /// @param cluster_front strip cluster on the first surface
  /// @param cluster_back strip cluster on the second surface
  /// @param theta the angle between the two strips
  /// @return global variances in R and Z directions.
  Acts::Vector2 getGlobalVars(const GeometryContext& gctx,
                              const cluster_t& cluster_front,
                              const cluster_t& cluster_back,
                              double theta) const;
};
}  // namespace Acts
#include "Acts/SpacePointFormation/detail/DoubleHitSpacePointBuilder.ipp"
