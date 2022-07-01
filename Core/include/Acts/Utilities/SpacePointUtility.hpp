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

/// @brief Storage container for variables related to the calculation of space
/// points
struct SpacePointParameters {
  /// Vector pointing from bottom to top end of first SDE
  Vector3 q;
  /// Vector pointing from bottom to top end of second SDE
  Vector3 r;
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3 s;
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3 t;
  /// Cross product between SpacePointParameters::q and
  /// SpacePointParameters::s
  Vector3 qs;
  /// Cross product between SpacePointParameters::r and
  /// SpacePointParameters::t
  Vector3 rt;
  /// Magnitude of SpacePointParameters::q
  double qmag = 0.;
  /// Parameter that determines the hit position on the first SDE
  double m = 0.;
  /// Parameter that determines the hit position on the second SDE
  double n = 0.;
  /// Regular limit of the absolut values of SpacePointParameters::m and
  /// SpacePointParameters::n
  double limit = 1.;
  /// Limit of SpacePointParameters::m and SpacePointParameters::n in case of
  /// variable vertex
  double limitExtended = 0.;
};

/// @class SpacePointUtility
///
class SpacePointUtility {
 public:
  using Measurement = Acts::BoundVariantMeasurement;

  /// Constructor
  SpacePointUtility(SpacePointBuilderConfig cfg) : m_config(cfg) {}

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
                          const Acts::Vector3& globalPos,
                          const Acts::SymMatrix2& localCov) const;

  /// @brief Calculate the global covariance from the front and back measurement in the strip SP formation
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param measFront The measurement on the front layer
  /// @param measBack The measurement on the back layer
  /// @param theta The angle between the two strips
  /// @return (rho, z) components of the global covariance
  Acts::Vector2 calcGlobalVars(const Acts::GeometryContext& gctx,
                               const Measurement& measFront,
                               const Measurement& measBack,
                               const Vector3& globalPos,
                               const double theta) const;

  bool calculateStripSPPosition(const std::pair<Vector3, Vector3>& stripEnds1,
                                const std::pair<Vector3, Vector3>& stripEnds2,
                                const Vector3& posVertex,
                                SpacePointParameters& spParams,
                                const double stripLengthTolerance) const;

  bool recoverSpacePoint(SpacePointParameters& spParams,
                         double stripLengthGapTolerance) const;

  double differenceOfMeasurementsChecked(const Vector3& pos1,
                                         const Vector3& pos2,
                                         const Vector3& posVertex,
                                         const double maxDistance,
                                         const double maxAngleTheta2,
                                         const double maxAnglePhi2) const;

  double calcPerpendicularProjection(
      const std::pair<Vector3, Vector3>& stripEnds1,
      const std::pair<Vector3, Vector3>& stripEnds2,
      SpacePointParameters& spParams) const;

 private:
  SpacePointBuilderConfig m_config;

  /// @brief Get the first component of the local covariance.
  /// @param meas The measurement
  /// @return the (0, 0) component of the local covariance
  double getLoc0Var(const Measurement& meas) const;
};

}  // namespace Acts
