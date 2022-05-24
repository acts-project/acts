// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"

#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

namespace Acts {

/// @class SpacePointUtility
///
class SpacePointUtility {
 public:
  using Measurement = Acts::BoundVariantMeasurement;

  /// Constructor
  SpacePointUtility(SpacePointBuilderConfig cfg) : m_config(cfg) {}

  /// @brief Getter method for the local coordinates of a measurement and its covariance
  ///
  /// @param meas measurement that holds the neccesary information of the hit position.
  /// @return vector of the local coordinates and covariance of the measurement on the surface
  std::pair<Acts::Vector2, Acts::SymMatrix2> getLocalPosCov(
      const Measurement& meas) const;

  /// @brief Getter method for the global coordinates of a measurement
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param meas measurement that holds the necessary
  /// information
  /// @return vectors of the global coordinates and covariance of the measurement
  std::pair<Vector3, Vector2> globalCoords(const GeometryContext& gctx,
                                           const Measurement& meas) const;

  /// @brief Get global covariance from the local position and covariance
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param geoId The geometry ID
  /// @param localPos The local position
  /// @param localCov The local covariance matrix
  /// @return (rho, z) components of the global covariance
  Acts::Vector2 globalCov(const Acts::GeometryContext& gctx,
                          const Acts::GeometryIdentifier& geoId,
                          const Acts::Vector2& localPos,
                          const Acts::SymMatrix2& localCov) const;

  /// @brief Get the first component of the local covariance.
  /// @param meas The measurement
  /// @return the (0, 0) component of the local covariance
  double getLoc0Var(const Measurement& meas) const;

  /// @brief Calculate the global covariance from the front and back measurement in the strip SP formation
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measFront The measurement on the front layer
  /// @param measBack The measurement on the back layer
  /// @param theta The angle between the two strips
  /// @return (rho, z) components of the global covariance
  Acts::Vector2 calcGlobalVars(const Acts::GeometryContext& gctx,
                               const Measurement& measFront,
                               const Measurement& measBack,
                               const double theta) const;

  /// @brief Get source link from the measurement
  /// @param meas The measurement
  const Acts::SourceLink* getSourceLink(const Measurement meas) const;

 private:
  SpacePointBuilderConfig m_config;
};

}  // namespace Acts
