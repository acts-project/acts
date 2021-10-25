// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @class OneHitSpacePointBuilder
///
/// After the particle interaction with surfaces are recorded and digitized
/// the clusters pixel detectors need further treatment. This class takes
/// the digitized clusters on a pixel detector element and provides the
/// corresponding space point.
///
template <typename spacepoint_t, typename cluster_t>
class SingleHitSpacePointBuilder {
 public:
  struct Config {
    // Tracking geometry for transformation lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  // Constructor
  SingleHitSpacePointBuilder(SingleHitSpacePointBuilderConfig cfg);
  ///// Default constructor
  SingleHitSpacePointBuilder() = default;

  /// @brief Calculates the space points out of a given collection of clusters
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clusters vector of clusters
  /// @param spacePointStorage storage of the results
  void calculateSpacePoints(const GeometryContext& gctx,
                            const std::vector<cluster_t>& clusters,
                            std::vector<spacepoint_t>& spacePointStorage) const;

 protected:
  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  ///
  /// @param clus cluster that holds the neccesary information of the 2D hit position.
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2 localCoords(const cluster_t& clus) const;

  /// @brief Getter method for the global coordinates of a cluster
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param clus cluster that holds the necessary
  /// information
  /// @return vector of the global coordinates and covariance of the cluster
  std::pair<Vector3, Vector2> globalCoords(const GeometryContext& gctx,
                                           const cluster_t& meas) const;

  // configuration of the single hit space point builder
  SingleHitSpacePointBuilderConfig m_config;
};
}  // namespace Acts
#include "Acts/SpacePointFormation/detail/SingleHitSpacePointBuilder.ipp"
