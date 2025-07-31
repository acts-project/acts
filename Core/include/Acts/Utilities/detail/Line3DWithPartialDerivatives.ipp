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
void Line3DWithPartialDerivatives<T>::updateParameters(const PaautoewPars) {
  constexpr auto x0 = static_cast<std::uint8_t>(ParIndex::x0);
  constexpr auto y0 = static_cast<std::uint8_t>(ParIndex::y0);
  constexpr auto theta = static_cast<std::uint8_t>(ParIndex::theta);
  constexpr auto phi = static_cast<std::uint8_t>(ParIndex::phi);

  m_pos[Acts::eX] = newPars[x0];
  m_pos[Acts::eY] = newPars[y0];

  const T cosTheta = std::cos(newPars[theta]);
  const T sinTheta = std::sin(newPars[theta]);
  const T cosPhi = std::cos(newPars[phi]);
  const T sinPhi = std::sin(newPars[phi]);

  m_dir = Vector{cosPhi * sinTheta, sinPhi * sinTheta, cosTheta};

  m_gradient[y0] = Vector::UnitY();
  m_gradient[x0] = Vector::UnitX();
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
  m_gradient[theta] = Vector{cosPhi * cosTheta, sinPhi * cosTheta, -sinTheta};
  m_gradient[phi] = Vector{-sinTheta * sinPhi, sinTheta * cosPhi, 0};
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
  constexpr auto idxThetaSq = vecIdxFromSymMat<s_nPars>(theta, theta);
  constexpr auto idxPhiSq = vecIdxFromSymMat<s_nPars>(phi, phi);
  constexpr auto idxPhiTheta = vecIdxFromSymMat<s_nPars>(theta, phi);
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
Line3DWithPartialDerivatives<T>::gradient(const ParIndex par) const {
  const std::uint8_t param = static_cast<std::uint8_t>(par);
  assert(param < m_gradient.size());
  return m_gradient[param];
}
template <std::floating_point T>
const Line3DWithPartialDerivatives<T>::Vector&
Line3DWithPartialDerivatives<T>::hessian(const ParIndex param1,
                                         const ParIndex param2) const {
  const auto idx{vecIdxFromSymMat<s_nPars>(static_cast<std::size_t>(param1),
                                           static_cast<std::size_t>(param2))};
  assert(idx < m_hessian.size());
  return m_hessian[idx];
}

template <std::floating_point T>
Line3DWithPartialDerivatives<T>::Vector Line3DWithPartialDerivatives<T>::point(
    const double lambda) const {
  return position() + lambda * direction();
}
}  // namespace Acts::detail
