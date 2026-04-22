// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"

namespace detray::detail {

/// @note Ported from Geant4 and simplified
template <concepts::scalar scalar_t>
struct density_effect_data {
  using scalar_type = scalar_t;

  constexpr density_effect_data() = default;

  DETRAY_HOST_DEVICE
  constexpr density_effect_data(const scalar_type a, const scalar_type m,
                                const scalar_type X0, const scalar_type X1,
                                const scalar_type I, const scalar_type nC,
                                const scalar_type delta0)
      : m_a(a),
        m_m(m),
        m_X0(X0),
        m_X1(X1),
        m_I(I * unit<scalar_type>::eV),
        m_nC(nC),
        m_delta0(delta0) {}

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  DETRAY_HOST_DEVICE constexpr bool operator==(
      const density_effect_data &rhs) const {
    return (m_a == rhs.get_A_density() && m_m == rhs.get_M_density() &&
            m_X0 == rhs.get_X0_density() && m_X1 == rhs.get_X1_density() &&
            m_I == rhs.get_mean_excitation_energy() &&
            m_nC == rhs.get_C_density() &&
            m_delta0 == rhs.get_delta0_density());
  }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_A_density() const { return m_a; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_M_density() const { return m_m; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_X0_density() const { return m_X0; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_X1_density() const { return m_X1; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_mean_excitation_energy() const { return m_I; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_C_density() const { return m_nC; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type get_delta0_density() const { return m_delta0; }

  /// Fitting parameters of Eq. 33.7 of RPP 2018
  /// @{
  scalar_type m_a = 0.f;
  scalar_type m_m = 0.f;
  scalar_type m_X0 = 0.f;
  scalar_type m_X1 = 0.f;
  /// @}
  /// Mean excitation energy in eV
  scalar_type m_I = 0.f;
  /// -C
  scalar_type m_nC = 0.f;
  /// Density-effect value delta(X_0)
  scalar_type m_delta0 = 0.f;
};

}  // namespace detray::detail
