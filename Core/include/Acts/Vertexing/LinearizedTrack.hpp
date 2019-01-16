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

class LinearizedTrack
{
public:
  /// Default constructor
  LinearizedTrack();

  /// Copy constructor
  LinearizedTrack(const LinearizedTrack& other) = default;

  /**
   * Constructor taking perigee parameters and covariance matrix
   * of track propagated to closest approach (PCA) of linearization point,
   * position and momentum Jacobian and const term.
   */
  LinearizedTrack(const Acts::ActsVectorD<5>&    paramsAtPCA,
                  const Acts::ActsSymMatrixD<5>& parCovarianceAtPCA,
                  const Acts::Vector3D&          linPoint,
                  const Acts::ActsMatrixD<5, 3>& positionJacobian,
                  const Acts::ActsMatrixD<5, 3>& momentumJacobian,
                  const Acts::Vector3D&       positionAtPCA,
                  const Acts::Vector3D&       momentumAtPCA,
                  const Acts::ActsVectorD<5>& constTerm);

  /// Assignment operator
  LinearizedTrack&
  operator=(const LinearizedTrack& other);

  /// Default constructor
  virtual ~LinearizedTrack() = default;

  const Acts::ActsVectorD<5>&
  parametersAtPCA() const
  {
    return m_paramsAtPCA;
  }

  const Acts::ActsSymMatrixD<5>&
  covarianceAtPCA() const
  {
    return m_parCovarianceAtPCA;
  }

  const Acts::Vector3D&
  linearizationPoint() const
  {
    return m_linPoint;
  }

  const Acts::ActsMatrixD<5, 3>&
  positionJacobian() const
  {
    return m_positionJacobian;
  }

  const Acts::ActsMatrixD<5, 3>&
  momentumJacobian() const
  {
    return m_momentumJacobian;
  }

  const Acts::Vector3D&
  positionAtPCA() const
  {
    return m_positionAtPCA;
  }

  const Acts::Vector3D&
  momentumAtPCA() const
  {
    return m_momentumAtPCA;
  }

  const Acts::ActsVectorD<5>&
  constantTerm() const
  {
    return m_constTerm;
  }

private:
  Acts::ActsVectorD<5>    m_paramsAtPCA;
  Acts::ActsSymMatrixD<5> m_parCovarianceAtPCA;
  Acts::Vector3D          m_linPoint;
  Acts::ActsMatrixD<5, 3> m_positionJacobian;
  Acts::ActsMatrixD<5, 3> m_momentumJacobian;
  Acts::Vector3D       m_positionAtPCA;
  Acts::Vector3D       m_momentumAtPCA;
  Acts::ActsVectorD<5> m_constTerm;
};

}  // namespace Acts
