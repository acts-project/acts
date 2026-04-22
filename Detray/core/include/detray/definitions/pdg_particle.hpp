// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"

// System include(s).
#include <cstdint>

namespace detray {

template <concepts::scalar scalar_t>
struct pdg_particle {
  using scalar_type = scalar_t;

  template <typename T>
  DETRAY_HOST_DEVICE constexpr pdg_particle(const std::int32_t pdg_num,
                                            const T mass, const T charge)
      : m_pdg_num(pdg_num),
        m_mass(static_cast<scalar_t>(mass)),
        m_charge(static_cast<scalar_t>(charge)) {}

  DETRAY_HOST_DEVICE
  constexpr std::int32_t pdg_num() const { return m_pdg_num; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type mass() const { return m_mass; }

  DETRAY_HOST_DEVICE
  constexpr scalar_type charge() const { return m_charge; }

  DETRAY_HOST_DEVICE
  constexpr bool is_valid() const { return m_pdg_num != 0; }

 private:
  std::int32_t m_pdg_num;
  scalar_type m_mass;
  scalar_type m_charge;
};

/// Apply the charge conjugation operator to a particle hypothesis @param ptc
template <concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr pdg_particle<scalar_t> charge_conjugation(
    const pdg_particle<scalar_t>& ptc) {
  return (ptc.charge() != 0)
             ? detray::pdg_particle<scalar_t>{-ptc.pdg_num(), ptc.mass(),
                                              -ptc.charge()}
             : ptc;
}

/// @returns an updated particle hypothesis according to the track qop
template <concepts::scalar scalar_t, typename track_t>
DETRAY_HOST_DEVICE constexpr pdg_particle<scalar_t> update_particle_hypothesis(
    const pdg_particle<scalar_t>& ptc, const track_t& params) {
  return (ptc.charge() * params.qop() > 0.f) ? ptc : charge_conjugation(ptc);
}

// Macro for declaring the particle
#define DETRAY_DECLARE_PARTICLE(PARTICLE_NAME, PDG_NUM, MASS, CHARGE) \
  template <concepts::scalar scalar_t>                                \
  struct PARTICLE_NAME final : public pdg_particle<scalar_t> {        \
    using base_type = pdg_particle<scalar_t>;                         \
    DETRAY_HOST_DEVICE                                                \
    constexpr PARTICLE_NAME() : base_type(PDG_NUM, MASS, CHARGE) {}   \
  }

// Declare PDG number 0 to be invalid
DETRAY_DECLARE_PARTICLE(invalid, 0, 0.f, 0.f);
// Declare some predefined particles
DETRAY_DECLARE_PARTICLE(electron, 11, constant<float>::m_e,
                        -1.f * unit<float>::e);
DETRAY_DECLARE_PARTICLE(positron, -11, constant<float>::m_e,
                        1.f * unit<float>::e);
// Muon mass from: https://physics.nist.gov/cgi-bin/cuu/Value?mmuc2mev (Visited
// Aug 1st, 2024)
DETRAY_DECLARE_PARTICLE(muon, 13, 105.6583755f * unit<float>::MeV,
                        -1.f * unit<float>::e);
DETRAY_DECLARE_PARTICLE(antimuon, -13, 105.6583755f * unit<float>::MeV,
                        1.f * unit<float>::e);
// Pion mass from: PDG 2024
DETRAY_DECLARE_PARTICLE(pion_zero, 111, 134.9768 * unit<float>::MeV, 0.f);
DETRAY_DECLARE_PARTICLE(pion_plus, 211, 139.57039f * unit<float>::MeV,
                        1.f * unit<float>::e);
DETRAY_DECLARE_PARTICLE(pion_minus, -211, 139.57039f * unit<float>::MeV,
                        -1.f * unit<float>::e);
DETRAY_DECLARE_PARTICLE(photon, 22, 0.f, 0.f);

}  // namespace detray
