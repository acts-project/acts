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

// Detray include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/pdg_particle.hpp"
#include "detray/definitions/units.hpp"
#include "detray/material/material.hpp"

// System include(s)
#include <cassert>

namespace detray::detail {

template <concepts::scalar scalar_t>
struct relativistic_quantities {
  using scalar_type = scalar_t;

  // values from RPP2018 table 33.1
  // Bethe formular prefactor. 1/mol unit is just a factor 1 here.
  static constexpr auto K{static_cast<scalar_type>(
      0.307075 * unit<double>::MeV * unit<double>::cm2)};
  // Energy scale for plasma energy.
  static constexpr auto PlasmaEnergyScale{
      static_cast<scalar_type>(28.816 * unit<double>::eV)};

  scalar_type m_qOverP{0.f};
  scalar_type m_q2OverBeta2{0.f};
  scalar_type m_beta{0.f};
  scalar_type m_beta2{0.f};
  scalar_type m_betaGamma{0.f};
  scalar_type m_gamma{0.f};
  scalar_type m_gamma2{0.f};
  scalar_type m_E{0.f};
  scalar_type m_Wmax{0.f};

  DETRAY_HOST_DEVICE
  relativistic_quantities(const pdg_particle<scalar_type>& ptc,
                          const scalar_type qop)
      : relativistic_quantities(ptc.mass(), qop, ptc.charge()) {};

  DETRAY_HOST_DEVICE relativistic_quantities(const scalar_type mass,
                                             const scalar_type qOverP,
                                             const scalar_type q)
      : m_qOverP{qOverP},
        // beta²/q² = (p/E)²/q² = p²/(q²m² + q²p²) = 1/(q² + (m²(q/p)²)
        // q²/beta² = q² + m²(q/p)²
        m_q2OverBeta2{q * q + (mass * qOverP) * (mass * qOverP)} {
    assert(m_qOverP != 0.f);
    assert(mass != 0.f);

    // 1/p = q/(qp) = (q/p)/q
    const scalar_type mOverP{
        mass * ((q != 0.f) ? math::fabs(qOverP / q) : math::fabs(qOverP))};
    const scalar_type pOverM{1.f / mOverP};
    // beta² = p²/E² = p²/(m² + p²) = 1/(1 + (m/p)²)
    m_beta2 = 1.f / (1.f + mOverP * mOverP);
    m_beta = math::sqrt(m_beta2);
    // beta*gamma = (p/sqrt(m² + p²))*(sqrt(m² + p²)/m) = p/m
    m_betaGamma = pOverM;
    // gamma = sqrt(m² + p²)/m = sqrt(1 + (p/m)²)
    m_gamma = math::sqrt(1.f + pOverM * pOverM);
    m_gamma2 = m_gamma * m_gamma;
    // E = gamma * mass;
    m_E = m_gamma * mass;

    const scalar_type mfrac{constant<scalar_type>::m_e / mass};

    // Wmax = 2m_e c^2 beta^2 gamma^2 / (1+2gamma*m_e/M + (m_e/M)^2)
    m_Wmax = (2.f * constant<scalar_type>::m_e * m_betaGamma * m_betaGamma) /
             (1.f + 2.f * m_gamma * mfrac + mfrac * mfrac);
  }

  /// @return 2 * mass * (beta * gamma)² mass term.
  DETRAY_HOST_DEVICE constexpr scalar_type compute_mass_term(
      const scalar_type mass) const {
    return 2.f * mass * m_betaGamma * m_betaGamma;
  }

  /// @return [(K/2) * (Z/A) * z^2 / beta^2 * density] in [energy/length]
  /// @brief defined in 34.12 of 2023 PDG review
  DETRAY_HOST_DEVICE constexpr scalar_type compute_epsilon_per_length(
      const material<scalar_type>& mat) const {
    return 0.5f * K * mat.molar_electron_density() * m_q2OverBeta2;
  }

  /// @return  (K/2) * (Z/A) * z^2 / beta^2 * density * path_length
  /// @brief defined in 34.12 of 2023 PDG review
  DETRAY_HOST_DEVICE constexpr scalar_type compute_epsilon(
      const material<scalar_type>& mat, const scalar_type path_length) const {
    return compute_epsilon_per_length(mat) * path_length;
  }

  /// @return the bethe_log_term of bethe equation,
  /// @brief bethe_log_term = 1/2 ln ( 2m_e c^2 beta^2 gamma^2 W_max/ I^2)
  DETRAY_HOST_DEVICE scalar_type
  compute_bethe_bloch_log_term(const material<scalar_type>& mat) const {
    const scalar_type I = mat.mean_excitation_energy();

    assert(I != 0.f);

    // u = 2 * m_e c^2* beta^2 * gamma^2
    const scalar_t u{compute_mass_term(constant<scalar_t>::m_e)};
    const scalar_type A = 0.5f * math::log(u * m_Wmax / (I * I));
    return A;
  }

  /// @return d(bethe_log_term)/dqop
  /// @brief  dA/dqop = - 1 / (2 * qop) * [4 - W_max/ (gamma M c^2) ]
  DETRAY_HOST_DEVICE scalar_type derive_bethe_bloch_log_term() const {
    assert(m_gamma != 0.f);
    assert(m_E != 0.f);
    const scalar_type dAdqop = -1.f / (2.f * m_qOverP) * (4.f - m_Wmax / m_E);
    return dAdqop;
  }

  /// @return d(beta^2)/dqop = - 2beta^2 / (qop * gamma^2)
  DETRAY_HOST_DEVICE scalar_type derive_beta2() const {
    assert(m_qOverP != 0.f);
    assert(m_gamma2 != 0.f);
    return -2.f * m_beta2 / (m_qOverP * m_gamma2);
  }

  /// @return the half of density correction factor (delta/2).
  DETRAY_HOST_DEVICE inline scalar_type compute_delta_half(
      const material<scalar_type>& mat) const {
    if (!mat.has_density_effect_data()) {
      // Apply the approximated density effector correction to tracks with
      // beta*gamma > 10
      //
      // @NOTE The approximated function is only accurate for beta*gamma >
      // 100 therefore, the cutoff with 10 might not be a good choice. The
      // cutoff also introduces a step where the energy loss function is
      // not smooth anymore.
      //
      // For further discussion, please follow the ATLAS JIRA Tickets:
      // ATLASRECTS-3144 and ATLASRECTS-7586 (ATLAS-restricted)
      if (m_betaGamma < 10.f) {
        return 0.f;
      }
      // Equation 34.6 of PDG2022
      // @NOTE A factor of 1000 is required to convert the unit of density
      // (mm^-3 to cm^-3)
      const scalar_type plasmaEnergy{
          PlasmaEnergyScale *
          math::sqrt(1000.f * mat.molar_electron_density())};
      return math::log(m_betaGamma * plasmaEnergy /
                       mat.mean_excitation_energy()) -
             0.5f;
    } else {
      const auto& density = mat.density_effect_data();

      const scalar_type cden{density.get_C_density()};
      const scalar_type mden{density.get_M_density()};
      const scalar_type aden{density.get_A_density()};
      const scalar_type x0den{density.get_X0_density()};
      const scalar_type x1den{density.get_X1_density()};

      const scalar_type x{math::log10(m_betaGamma)};

      scalar_type delta{0.f};

      // From Geant4
      // processes/electromagnetic/lowenergy/src/G4hBetheBlochModel.cc
      if (x < x0den) {
        delta = 0.f;
        // @TODO: Add a branch for conductors (Eq 34.7 of
        // https://pdg.lbl.gov/2023/reviews/rpp2023-rev-particle-detectors-accel.pdf)
      } else {
        delta = 2.f * constant<scalar_type>::ln10 * x - cden;
        if (x < x1den)
          delta += aden * math::pow((x1den - x), mden);
      }

      return 0.5f * delta;
    }
  }

  /// @return the derivation of the density correction factor delta/2.
  DETRAY_HOST_DEVICE inline scalar_type derive_delta_half(
      const material<scalar_type>& mat) const {
    assert(m_qOverP != 0.f);

    if (!mat.has_density_effect_data()) {
      // d(ln(betagamma))/dqop = -1/qop
      return -1.f / m_qOverP;
    } else {
      const auto& density = mat.density_effect_data();

      // const scalar_type cden{density.get_C_density()};
      const scalar_type mden{density.get_M_density()};
      const scalar_type aden{density.get_A_density()};
      const scalar_type x0den{density.get_X0_density()};
      const scalar_type x1den{density.get_X1_density()};

      const scalar_type x{math::log10(m_betaGamma)};

      scalar_type delta{0.f};

      // From Geant4
      // processes/electromagnetic/lowenergy/src/G4hBetheBlochModel.cc
      if (x < x0den) {
        delta = 0.f;
        // @TODO: Add a branch for conductors (Eq 34.7 of
        // https://pdg.lbl.gov/2023/reviews/rpp2023-rev-particle-detectors-accel.pdf)
      } else {
        delta = -2.f / m_qOverP;
        if (x < x1den) {
          delta += aden * mden / (m_qOverP * constant<scalar_type>::ln10) *
                   math::pow(x1den - x, mden - 1.f);
        }
      }

      return 0.5f * delta;
    }
  }
};

}  // namespace detray::detail
