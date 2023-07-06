// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/Charge.hpp"

#include <utility>

namespace Acts {

/// @brief Particle hypothesis used in reconstruction
///
/// Our reconstruction hypothesis consists of absolute PDG code, mass and
/// absolute charge
template <typename charge_t>
class GenericParticleHypothesis {
 public:
  using ChargeType = charge_t;

  constexpr GenericParticleHypothesis(PdgParticle absPdg, float mass,
                                      ChargeType chargeType)
      : m_absPdg{absPdg}, m_mass{mass}, m_chargeType{std::move(chargeType)} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  template <typename other_charge_t>
  constexpr GenericParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : m_absPdg{other.absPdg()},
        m_mass{other.mass()},
        m_chargeType{other.chargeType()} {}

  constexpr PdgParticle absPdg() const noexcept { return m_absPdg; }

  constexpr float mass() const noexcept { return m_mass; }

  constexpr float absoluteCharge() const noexcept {
    return m_chargeType.absQ();
  }

  template <typename T>
  constexpr auto extractCharge(T qOverP) const noexcept {
    return m_chargeType.extractCharge(qOverP);
  }

  template <typename T>
  constexpr auto extractMomentum(T qOverP) const noexcept {
    return m_chargeType.extractMomentum(qOverP);
  }

  template <typename P, typename Q>
  constexpr auto chargeOverMomentum(P p, Q q) const noexcept {
    return m_chargeType.chargeOverMomentum(p, q);
  }

  constexpr const ChargeType& chargeType() const noexcept {
    return m_chargeType;
  }

 private:
  PdgParticle m_absPdg;
  float m_mass;
  ChargeType m_chargeType;

  friend bool operator==(const GenericParticleHypothesis<ChargeType>& lhs,
                         const GenericParticleHypothesis<ChargeType>& rhs) {
    return (lhs.m_absPdg == rhs.m_absPdg) && (lhs.m_mass == rhs.m_mass) &&
           (lhs.m_chargeType == rhs.m_chargeType);
  }
  friend bool operator!=(const GenericParticleHypothesis<ChargeType>& lhs,
                         const GenericParticleHypothesis<ChargeType>& rhs) {
    return (lhs.m_absPdg != rhs.m_absPdg) || (lhs.m_mass != rhs.m_mass) ||
           (lhs.m_chargeType != rhs.m_chargeType);
  }
};

}  // namespace Acts
