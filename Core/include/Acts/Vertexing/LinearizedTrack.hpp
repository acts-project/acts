// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class LinearizedTrack
///
/// Class for linear expansion of track parameters in vicinity of vertex
///
/// The measurement equation is linearized in the following way:
///
/// F_k= D_k (x_k - x_0k) + E_k (p_k - p_0k) + F^0_k
///
/// where F_k are the parameters at perigee nearest to the linearization point,
/// x_k is the position of the vertex, p_k the track momentum at the vertex,
/// and F^0_k is the constant term of expansion. D_k and E_k are matrices
/// of derivatives, denoted hereafter as "positionJacobian" and
/// "momentumJacobian" respectively.
///

struct LinearizedTrack {
  /// @brief Constructor taking perigee parameters and covariance matrix
  /// of track propagated to closest approach (PCA) of linearization point,
  /// position and momentum Jacobian and const term.
  ///
  /// @param paramsAtPCA Parameters at point of closest approach
  /// @param parCovarianceAtPCA Parameter covariance matrix at point of closest
  /// approach
  /// @param linPoint Linearization point
  /// @param posJacobian Position jacobian
  /// @param momJacobian Momentum jacobian
  /// @param position Position at point of closest approach
  /// @param momentum Momentum at point of closest approach
  /// @param constTerm Constant term in taylor expansion

  LinearizedTrack() = default;

  LinearizedTrack(const TrackVector& paramsAtPCA,
                  const TrackSymMatrix& parCovarianceAtPCA,
                  const Vector3D& linPoint,
                  const SpacePointToTrackMatrix& posJacobian,
                  const SpacePointToTrackMatrix& momJacobian,
                  const Vector3D& position, const Vector3D& momentum,
                  const TrackVector& constTerm)
      : parametersAtPCA(paramsAtPCA),
        covarianceAtPCA(parCovarianceAtPCA),
        linearizationPoint(linPoint),
        positionJacobian(posJacobian),
        momentumJacobian(momJacobian),
        positionAtPCA(position),
        momentumAtPCA(momentum),
        constantTerm(constTerm) {}

  TrackVector parametersAtPCA{TrackVector::Zero()};
  TrackSymMatrix covarianceAtPCA{TrackSymMatrix::Zero()};
  Vector3D linearizationPoint{Vector3D::Zero()};
  SpacePointToTrackMatrix positionJacobian{SpacePointToTrackMatrix::Zero()};
  SpacePointToTrackMatrix momentumJacobian{SpacePointToTrackMatrix::Zero()};
  Vector3D positionAtPCA{Vector3D::Zero()};
  Vector3D momentumAtPCA{Vector3D::Zero()};
  TrackVector constantTerm{TrackVector::Zero()};
};

}  // namespace Acts
