// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/**
 * @class LinearizedTrack.h
 *
 * Class for linear expansion of track parameters in vicinity of vertex
 *
 * The measurement equation is linearized in the following way:
 *
 * F_k= D_k (x_k - x_0k) + E_k (p_k - p_0k) + F^0_k
 *
 * where F_k are the parameters at perigee nearest to the lin point,
 * x_k is the position of the vertex, p_k the track momentum at the vertex,
 * and F^0_k is the constant term of expansion. D_k and E_k are matrices
 * of derivatives, denoted hereafter as "positionJacobian" and
 * "momentumJacobian"
 * respectively.
 */

struct LinearizedTrack
{

  /**
   * Constructor taking perigee parameters and covariance matrix
   * of track propagated to closest approach (PCA) of linearization point,
   * position and momentum Jacobian and const term.
   */
  LinearizedTrack(const ActsVectorD<5>&    paramsAtPCA,
                  const ActsSymMatrixD<5>& parCovarianceAtPCA,
                  const Vector3D&          linPoint,
                  const ActsMatrixD<5, 3>& positionJacobian,
                  const ActsMatrixD<5, 3>& momentumJacobian,
                  const Vector3D&       positionAtPCA,
                  const Vector3D&       momentumAtPCA,
                  const ActsVectorD<5>& constTerm);


  Acts::ActsVectorD<5>    parametersAtPCA;
  Acts::ActsSymMatrixD<5> covarianceAtPCA;
  Acts::Vector3D          linearizationPoint;
  Acts::ActsMatrixD<5, 3> positionJacobian;
  Acts::ActsMatrixD<5, 3> momentumJacobian;
  Acts::Vector3D       positionAtPCA;
  Acts::Vector3D       momentumAtPCA;
  Acts::ActsVectorD<5> constantTerm;
};

}  // namespace Acts
