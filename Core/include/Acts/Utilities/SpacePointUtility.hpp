// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.hpp"
#include "Acts/Utilities/Result.hpp"

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
  Vector3 firstBtmToTop;
  /// Vector pointing from bottom to top end of second SDE
  Vector3 secondBtmToTop;
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3 vtxToFirstMid2;
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3 vtxToSecondMid2;
  /// Cross product between firstBtmToTop and vtxToFirstMid2
  Vector3 firstBtmToTopXvtxToFirstMid2;
  /// Cross product between secondBtmToTop and vtxToSecondMid2
  Vector3 secondBtmToTopXvtxToSecondMid2;
  /// Magnitude of SpacePointParameters::firstBtmToTop
  double mag_firstBtmToTop = 0.;
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
  /// Constructor
  SpacePointUtility(SpacePointBuilderConfig cfg) : m_config(std::move(cfg)) {}

  /// @brief Getter method for the global coordinates of a SourceLink
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param slink SourceLink that holds the necessary
  /// information
  /// @param par local position
  /// @param cov local covariance
  /// @return vectors of the global coordinates and covariance of the SourceLink
  std::pair<Vector3, Vector2> globalCoords(const GeometryContext& gctx,
                                           const SourceLink& slink,
                                           const BoundVector& par,
                                           const BoundSymMatrix& cov) const;

  /// @brief Get rho and z covariance from the local position and covariance
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param geoId The geometry ID
  /// @param globalPos The global position
  /// @param localCov The local covariance matrix
  /// @return (rho, z) components of the global covariance

  Vector2 rhoZCovariance(const GeometryContext& gctx,
                         const GeometryIdentifier& geoId,
                         const Vector3& globalPos,
                         const SymMatrix2& localCov) const;

  /// @brief Calculate the rho and z covariance from the front and back SourceLink in the strip SP formation
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param slinkFront The SourceLink on the front layer
  /// @param slinkBack The SourceLink on the back layer
  /// @param paramCovAccessor function to extract local position and covariance from SourceLink
  /// @param globalPos global position
  /// @param theta The angle between the two strips
  /// @return (rho, z) components of the global covariance
  Vector2 calcRhoZVars(
      const GeometryContext& gctx, const SourceLink& slinkFront,
      const SourceLink& slinkBack,
      const std::function<std::pair<const BoundVector, const BoundSymMatrix>(
          SourceLink)>& paramCovAccessor,
      const Vector3& globalPos, const double theta) const;

  /// @brief This function performs a straight forward calculation of a space
  /// point and returns whether it was succesful or not.
  ///
  /// @param [in] stripEnds1 Top and bottom end of the first strip
  /// @param [in] stripEnds2 Top and bottom end of the second strip
  /// @param [in] posVertex Position of the vertex
  /// @param [in, out] spParams Data container of the calculations
  /// @param [in] stripLengthTolerance Tolerance scaling factor on the strip
  /// detector element length
  ///
  /// @return Result whether the space point calculation was succesful
  Result<void> calculateStripSPPosition(
      const std::pair<Vector3, Vector3>& stripEnds1,
      const std::pair<Vector3, Vector3>& stripEnds2, const Vector3& posVertex,
      SpacePointParameters& spParams, const double stripLengthTolerance) const;

  /// @brief This function tests if a space point can be estimated by a more
  /// tolerant treatment of construction. In fact, this function indirectly
  /// allows shifts of the vertex.
  ///
  /// @param [in] spParams container that stores geometric parameters and rules of
  /// the space point formation
  /// @param [in] stripLengthGapTolerance Tolerance scaling factor of the gap
  /// between strip detector elements
  ///
  /// @return indicator if the test was successful
  Result<void> recoverSpacePoint(SpacePointParameters& spParams,
                                 double stripLengthGapTolerance) const;

  /// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two SourceLinks
  ///
  /// @param [in] pos1 position of the first SourceLink
  /// @param [in] pos2 position the second SourceLink
  /// @param [in] posVertex Position of the vertex
  /// @param [in] maxDistance Maximum distance between two SourceLinks
  /// @param [in] maxAngleTheta2 Maximum squared theta angle between two
  /// SourceLinks
  /// @param [in] maxAnglePhi2 Maximum squared phi angle between two SourceLinks
  ///
  /// @return Result with the squared sum within configuration parameters.
  Result<double> differenceOfMeasurementsChecked(
      const Vector3& pos1, const Vector3& pos2, const Vector3& posVertex,
      const double maxDistance, const double maxAngleTheta2,
      const double maxAnglePhi2) const;

  /// @brief Calculates a space point whithout using the vertex
  /// @note This is mostly to resolve space points from cosmic data
  /// @param stripEnds1 The ends of one strip
  /// @param stripEnds2 The ends of another strip
  /// @param spParams SpacePointParamaters for the SP
  /// @return parameter that indicates the location of the space point; returns
  /// 1. if it failed
  /// @note The meaning of the parameter is explained in more detail in the
  /// function body
  Result<double> calcPerpendicularProjection(
      const std::pair<Vector3, Vector3>& stripEnds1,
      const std::pair<Vector3, Vector3>& stripEnds2,
      SpacePointParameters& spParams) const;

 private:
  SpacePointBuilderConfig m_config;
  std::error_code m_error;
};

}  // namespace Acts
