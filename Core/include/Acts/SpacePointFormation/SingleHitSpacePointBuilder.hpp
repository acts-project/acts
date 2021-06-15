// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
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
template <typename spacepoint_t, typename source_link_t>
class SingleHitSpacePointBuilder {
 public:
  struct Config {
    // Tracking geometry for transformation lookup.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };

  // Constructor
  SingleHitSpacePointBuilder(SpacePointBuilderConfig cfg);
  ///// Default constructor
  SingleHitSpacePointBuilder() = default;

  /// @brief Calculates the space points out of a given collection of clusters
  /// and stores the results
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measurements vector of measurements
  /// @param spacePointStorage storage of the results
  void calculateSpacePoints(
      const GeometryContext& gctx,
      const std::vector<BoundVariantMeasurement<source_link_t>>& measurements,
      std::vector<spacepoint_t>& spacePointStorage) const;

 protected:
  /// @brief Getter method for the local coordinates of a cluster
  /// on its corresponding surface
  ///
  /// @param meas object related to the measurement that holds the necessary
  /// information
  /// @return vector of the local coordinates of the cluster on the surface
  Vector2 localCoords(const BoundVariantMeasurement<source_link_t>& meas) const;

  /// @brief Getter method for the global coordinates of a cluster
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param meas object related to the measurement that holds the necessary
  /// information
  /// @return vector of the global coordinates of the cluster
  std::pair<Vector3, Vector2> globalCoords(
      const GeometryContext& gctx,
      const BoundVariantMeasurement<source_link_t>& meas) const;

  SpacePointBuilderConfig m_config;
};
}  // namespace Acts
#include "Acts/SpacePointFormation/detail/SingleHitSpacePointBuilder.ipp"
