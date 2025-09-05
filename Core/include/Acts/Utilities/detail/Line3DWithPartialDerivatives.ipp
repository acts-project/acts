// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace Acts::detail {

template <std::floating_point T>
Line3DWithPartialDerivatives<T>::ParamVector
Line3DWithPartialDerivatives<T>::parameters() const {
  ParamVector pars{};
  pars[toUnderlying(ParIndex::x0)] = m_pos[Acts::eX];
  pars[toUnderlying(ParIndex::y0)] = m_pos[Acts::eY];
  pars[toUnderlying(ParIndex::theta)] = Acts::VectorHelpers::theta(m_dir);
  pars[toUnderlying(ParIndex::phi)] = Acts::VectorHelpers::phi(m_dir);
  return pars;
}

template <std::floating_point T>
template <std::size_t N>
Line3DWithPartialDerivatives<T>::Line3DWithPartialDerivatives(
    const std::array<T, N>& initPars) noexcept {
  updateParameters(initPars);
}

template <std::floating_point T>
template <std::size_t N>
void Line3DWithPartialDerivatives<T>::updateParameters(
    const std::array<T, N>& newPars) noexcept
  requires(N >= s_nPars)
{
  m_pos[Acts::eX] = newPars[toUnderlying(ParIndex::x0)];
  m_pos[Acts::eY] = newPars[toUnderlying(ParIndex::y0)];

  const T cosTheta = std::cos(newPars[toUnderlying(ParIndex::theta)]);
  const T sinTheta = std::sin(newPars[toUnderlying(ParIndex::theta)]);
  const T cosPhi = std::cos(newPars[toUnderlying(ParIndex::phi)]);
  const T sinPhi = std::sin(newPars[toUnderlying(ParIndex::phi)]);

  m_dir = Vector{cosPhi * sinTheta, sinPhi * sinTheta, cosTheta};

  m_gradient[toUnderlying(ParIndex::y0)] = Vector::UnitY();
  m_gradient[toUnderlying(ParIndex::x0)] = Vector::UnitX();
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
  m_gradient[toUnderlying(ParIndex::theta)] =
      Vector{cosPhi * cosTheta, sinPhi * cosTheta, -sinTheta};
  m_gradient[toUnderlying(ParIndex::phi)] =
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
  constexpr auto idxThetaSq = vecIdxFromSymMat<s_nPars>(
      toUnderlying(ParIndex::theta), toUnderlying(ParIndex::theta));
  constexpr auto idxPhiSq = vecIdxFromSymMat<s_nPars>(
      toUnderlying(ParIndex::phi), toUnderlying(ParIndex::phi));
  constexpr auto idxPhiTheta = vecIdxFromSymMat<s_nPars>(
      toUnderlying(ParIndex::theta), toUnderlying(ParIndex::phi));
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
  const auto param = static_cast<std::uint8_t>(par);
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
