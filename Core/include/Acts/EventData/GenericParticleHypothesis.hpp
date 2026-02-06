// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/EventData/ChargeConcept.hpp"

#include <cassert>
#include <ostream>
#include <utility>

namespace Acts {

/// @brief Particle hypothesis used in reconstruction
///
/// The reconstruction hypothesis consists of absolute PDG code, mass and
/// absolute charge.
template <ChargeConcept charge_t>
class GenericParticleHypothesis {
 public:
  /// Type alias for charge type used in particle hypothesis
  using ChargeType = charge_t;

  /// Creates a particle hypothesis using absolute PDG, mass and the charge
  /// type.
  ///
  /// @param absPdg the absolute PDG
  /// @param mass the particle mass
  /// @param chargeType the type of charge
  constexpr GenericParticleHypothesis(PdgParticle absPdg, float mass,
                                      ChargeType chargeType)
      : m_absPdg{absPdg}, m_mass{mass}, m_chargeType{std::move(chargeType)} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Creates a particle hypothesis using the absolute PDG.
  /// The mass and charge is looked up using @ref findMass and @ref findCharge.
  /// If the lookup fails an exception is thrown.
  ///
  /// @param absPdg the absolute PDG
  explicit GenericParticleHypothesis(PdgParticle absPdg)
      : m_absPdg{absPdg},
        m_mass{findMass(absPdg).value()},
        m_chargeType{std::abs(findCharge(absPdg).value())} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Copy from another charge hypothesis.
  /// @param other The other particle hypothesis to copy from
  template <typename other_charge_t>
  explicit constexpr GenericParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : m_absPdg{other.absolutePdg()},
        m_mass{other.mass()},
        m_chargeType{other.chargeType()} {}

  /// Get the hypothesized absolute PDG.
  /// @return The absolute PDG particle identifier
  constexpr PdgParticle absolutePdg() const noexcept { return m_absPdg; }

  /// Get the hypothesized mass.
  /// @return The particle mass in natural units
  constexpr float mass() const noexcept { return m_mass; }

  /// Get the hypothesized absolute charge.
  /// @return The absolute charge magnitude
  constexpr float absoluteCharge() const noexcept {
    return m_chargeType.absQ();
  }

  /// Extracts the signed charge from the `q over p` track parameter using the
  /// charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted signed charge
  template <typename T>
  constexpr auto extractCharge(T qOverP) const noexcept {
    return m_chargeType.extractCharge(qOverP);
  }

  /// Extracts the particle momentum from the `q over p` track parameter using
  /// the charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted absolute momentum
  template <typename T>
  constexpr auto extractMomentum(T qOverP) const noexcept {
    return m_chargeType.extractMomentum(qOverP);
  }

  /// Calculate the `q over p` track parameter with the given absolute momentum
  /// and charge.
  ///
  /// @param momentum the absolute momentum.
  /// @param signedQ the signed charge.
  /// @return The calculated charge over momentum ratio
  template <typename P, typename Q>
  constexpr auto qOverP(P momentum, Q signedQ) const noexcept {
    return m_chargeType.qOverP(momentum, signedQ);
  }

  /// Get the hypothesized charge type.
  /// @return Reference to the charge type object
  constexpr const ChargeType& chargeType() const noexcept {
    return m_chargeType;
  }

  /// Output stream representation of the particle hypothesis
  /// @param os Output stream to write to
  /// @return Modified output stream for chaining\n
  std::ostream& toStream(std::ostream& os) const {
    os << "ParticleHypothesis{absPdg=";
    if (auto shortString = pdgToShortAbsString(absolutePdg())) {
      os << *shortString;
    } else {
      os << absolutePdg();
    }
    os << ", mass=" << mass() << ", absCharge=" << absoluteCharge() << "}";
    return os;
  }

  /// Output stream operator for particle hypothesis
  /// @param os Output stream to write to
  /// @param particleHypothesis The particle hypothesis to output
  /// @return Reference to output stream for chaining
  friend std::ostream& operator<<(
      std::ostream& os, const GenericParticleHypothesis& particleHypothesis) {
    return particleHypothesis.toStream(os);
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
};

}  // namespace Acts
