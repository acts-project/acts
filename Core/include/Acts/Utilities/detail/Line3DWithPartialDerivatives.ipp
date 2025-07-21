// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

namespace Acts::detail {

template <std::floating_point T>
void Line3DWithPartialDerivatives<T>::updateParameters(
    const ParamVector& newPars) {
  m_pos[Acts::eX] = newPars[ParIndices::x0];
  m_pos[Acts::eY] = newPars[ParIndices::y0];

  const T cosTheta = std::cos(newPars[ParIndices::theta]);
  const T sinTheta = std::sin(newPars[ParIndices::theta]);
  const T cosPhi = std::cos(newPars[ParIndices::phi]);
  const T sinPhi = std::sin(newPars[ParIndices::phi]);

  m_dir = Vector{cosPhi * sinTheta, sinPhi * sinTheta, cosTheta};

  m_gradient[ParIndices::y0] = Vector::UnitY();
  m_gradient[ParIndices::x0] = Vector::UnitX();
  /**            x_{0}                cos (phi) sin (theta)
   *  position = y_{0}  , Direction = sin (phi) sin (theta)
   *             0                         cos theta
   *
   *
   *   d Direction     cos (phi) cos(theta)     d Direction     -sin(phi) sin
   *(theta)
   *  -------------=   sin (phi) cos(theta)     ------------ =   cos(phi) sin
   *(theta) dTheta         - sin (theta)             dPhi                 0
   *
   *******************************************************************************/
  m_gradient[ParIndices::theta] =
      Vector{cosPhi * cosTheta, sinPhi * cosTheta, -sinTheta};
  m_gradient[ParIndices::phi] =
      Vector{-sinTheta * sinPhi, sinTheta * cosPhi, 0};
  /*********************************************************************************
   *   Non-vanishing second order derivatives
   *
   *    d^{2} Direction                 d^{2} Direction                     cos
   *phi
   *    ------------- = - Direction ,   ------------      = - sin(theta)    sin
   *phi d^{2} theta                      d^{2} phi                             0
   *
   *   d^{2} Direction                 -sin phi
   *   -------------     = cos(theta)   cos phi
   *    d theta dPhi                       0
   ************************************************************************************/
  constexpr std::size_t idxThetaSq =
      vecIdxFromSymMat<ParIndices::nPars>(ParIndices::theta, ParIndices::theta);
  constexpr std::size_t idxPhiSq =
      vecIdxFromSymMat<ParIndices::nPars>(ParIndices::phi, ParIndices::phi);
  constexpr std::size_t idxPhiTheta =
      vecIdxFromSymMat<ParIndices::nPars>(ParIndices::theta, ParIndices::phi);
  m_hessian[idxThetaSq] = -m_dir;
  m_hessian[idxPhiSq] = -sinTheta * Vector{cosPhi, sinPhi, 0.};
  m_hessian[idxPhiTheta] = cosTheta * Vector{-sinPhi, cosPhi, 0.};
}

template <std::floating_point T>
const Line3DWithPartialDerivatives<T>::Vector&
Line3DWithPartialDerivatives<T>::position() const {
  return m_pos;
}
template <std::floating_point T>
const Line3DWithPartialDerivatives<T>::Vector&
Line3DWithPartialDerivatives<T>::direction() const {
  return m_dir;
}
template <std::floating_point T>
const Line3DWithPartialDerivatives<T>::Vector&
Line3DWithPartialDerivatives<T>::gradient(const std::size_t param) const {
  assert(param < m_gradient.size());
  return m_gradient[param];
}
template <std::floating_point T>
const Line3DWithPartialDerivatives<T>::Vector&
Line3DWithPartialDerivatives<T>::hessian(const std::size_t param1,
                                         const std::size_t param2) const {
  const std::size_t idx{vecIdxFromSymMat<s_nPars>(param1, param2)};
  assert(idx < m_hessian.size());
  return m_hessian[idx];
}

template <std::floating_point T>
Line3DWithPartialDerivatives<T>::Vector Line3DWithPartialDerivatives<T>::point(
    const double lambda) const {
  return position() + lambda * direction();
}
}  // namespace Acts::detail
