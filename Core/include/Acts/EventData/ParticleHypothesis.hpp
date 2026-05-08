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
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ChargeHypothesis.hpp"

#include <cassert>
#include <iostream>

namespace Acts {

// TODO In principle the factory methods could provide a reference to a static
// instance which would avoid copying the particle hypothesis and potentially
// save some memory. But constexpr+static seems to require C++2b extension.

/// @ingroup eventdata
/// @defgroup eventdata-particlehypothesis Particle hypothesis for track reconstruction
/// @{

/// @brief Particle hypothesis used in reconstruction
///
/// The reconstruction hypothesis consists of absolute PDG code, mass and
/// absolute charge.
class ParticleHypothesis final {
 public:
  /// Create a muon particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Muon particle hypothesis with any charge type
  static ParticleHypothesis muon(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::eMuon);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create a charged pion particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Charged pion particle hypothesis with any charge type
  static ParticleHypothesis pion(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::ePionPlus);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create an electron particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Electron particle hypothesis with any charge type
  static ParticleHypothesis electron(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::eElectron);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create a charged kaon particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Charged kaon particle hypothesis with any charge type
  static ParticleHypothesis kaon(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::eKaonPlus);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create a proton particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Proton particle hypothesis with any charge type
  static ParticleHypothesis proton(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::eProton);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }

  /// Create a photon particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Photon particle hypothesis with any charge type
  static ParticleHypothesis photon(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::eGamma);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create a neutral pion particle hypothesis
  /// @param momentum The optional particle momentum
  /// @return Neutral pion particle hypothesis with any charge type
  static ParticleHypothesis pion0(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::ePionZero);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }

  /// Create a pion-like particle hypothesis with custom charge
  /// @param absoluteCharge The absolute charge value
  /// @param momentum The optional particle momentum
  /// @return Pion-like particle hypothesis with any charge type
  static ParticleHypothesis pionLike(
      float absoluteCharge, std::optional<double> momentum = std::nullopt) {
    return ParticleHypothesis(pion().absolutePdg(), pion().mass(),
                              ChargeHypothesis{absoluteCharge}, momentum);
  }

  /// Create a neutral geantino particle hypothesis (massless neutral particle)
  /// @param momentum The optional particle momentum
  /// @return Neutral geantino particle hypothesis with any charge type
  static ParticleHypothesis geantino(
      std::optional<double> momentum = std::nullopt) {
    static const ParticleHypothesis cache(PdgParticle::eInvalid, 0,
                                          ChargeHypothesis{0});
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create a charged geantino particle hypothesis with unit charge
  /// @param momentum The optional particle momentum
  /// @return Charged geantino particle hypothesis with any charge type
  static ParticleHypothesis chargedGeantino(
      std::optional<double> momentum = std::nullopt) {
    static const auto cache = chargedGeantino(Acts::UnitConstants::e);
    ParticleHypothesis result = cache;
    result.m_momentum = momentum;
    return result;
  }
  /// Create a charged geantino particle hypothesis with custom charge
  /// @param absoluteCharge The absolute charge value
  /// @return Charged geantino particle hypothesis with any charge type
  static ParticleHypothesis chargedGeantino(
      float absoluteCharge, std::optional<double> momentum = std::nullopt) {
    ParticleHypothesis result = ParticleHypothesis(
        PdgParticle::eInvalid, 0, ChargeHypothesis{absoluteCharge});
    result.m_momentum = momentum;
    return result;
  }

  /// Creates a particle hypothesis using absolute PDG, mass and the charge
  /// type.
  ///
  /// @param absPdg the absolute PDG
  /// @param mass the particle mass
  /// @param absCharge the absolute charge
  /// @param momentum the optional particle momentum
  constexpr ParticleHypothesis(PdgParticle absPdg, float mass, float absCharge,
                               std::optional<double> momentum = std::nullopt)
      : m_absPdg{absPdg},
        m_mass{mass},
        m_charge{absCharge},
        m_momentum{momentum} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Creates a particle hypothesis using absolute PDG, mass and the charge
  /// type.
  ///
  /// @param absPdg the absolute PDG
  /// @param mass the particle mass
  /// @param charge the charge type
  /// @param momentum the optional particle momentum
  constexpr ParticleHypothesis(PdgParticle absPdg, float mass,
                               ChargeHypothesis charge,
                               std::optional<double> momentum = std::nullopt)
      : m_absPdg{absPdg}, m_mass{mass}, m_charge{charge}, m_momentum{momentum} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Creates a particle hypothesis using the absolute PDG.
  /// The mass and charge is looked up using @ref findMass and @ref findCharge.
  /// If the lookup fails an exception is thrown.
  ///
  /// @param absPdg the absolute PDG
  /// @param momentum the optional particle momentum
  explicit ParticleHypothesis(PdgParticle absPdg,
                              std::optional<double> momentum = std::nullopt)
      : m_absPdg{absPdg},
        m_mass{findMass(absPdg).value()},
        m_charge{std::abs(findCharge(absPdg).value())},
        m_momentum{momentum} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Get the hypothesized absolute PDG.
  /// @return The absolute PDG particle identifier
  constexpr PdgParticle absolutePdg() const noexcept { return m_absPdg; }

  /// Get the hypothesized mass.
  /// @return The particle mass
  constexpr float mass() const noexcept { return m_mass; }

  /// Get the hypothesized absolute charge.
  /// @return The absolute charge magnitude
  float absoluteCharge() const noexcept { return m_charge.absoluteCharge(); }

  /// Get the optional hypothesized momentum.
  /// @return The optional particle momentum
  std::optional<double> momentum() const noexcept { return m_momentum; }

  /// Extracts the signed charge from the `q over p` track parameter using the
  /// charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted signed charge
  constexpr float extractCharge(double qOverP) const noexcept {
    return m_charge.extractCharge(qOverP);
  }

  /// Extracts the particle momentum from the `q over p` track parameter using
  /// the charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted absolute momentum
  constexpr double extractMomentum(double qOverP) const noexcept {
    if (m_momentum.has_value()) {
      return *m_momentum;
    }
    return m_charge.extractMomentum(qOverP);
  }

  /// Calculate the `q over p` track parameter with the given absolute momentum
  /// and charge.
  ///
  /// @param momentum the absolute momentum.
  /// @param signedQ the signed charge.
  /// @return The calculated charge over momentum ratio
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    return m_charge.qOverP(momentum, signedQ);
  }

  /// Get the hypothesized charge.
  /// @return Reference to the charge type object
  constexpr const ChargeHypothesis& charge() const noexcept { return m_charge; }

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
    os << ", mass=" << mass() << ", absCharge=" << absoluteCharge();
    if (momentum().has_value()) {
      os << ", momentum=" << *momentum();
    }
    os << "}";
    return os;
  }

  /// Output stream operator for particle hypothesis
  /// @param os Output stream to write to
  /// @param particleHypothesis The particle hypothesis to output
  /// @return Reference to output stream for chaining
  friend std::ostream& operator<<(
      std::ostream& os, const ParticleHypothesis& particleHypothesis) {
    return particleHypothesis.toStream(os);
  }

 private:
  PdgParticle m_absPdg{PdgParticle::eInvalid};
  float m_mass{0};
  ChargeHypothesis m_charge;
  std::optional<double> m_momentum;

  friend bool operator==(const ParticleHypothesis& lhs,
                         const ParticleHypothesis& rhs) = default;
};

/// @}

}  // namespace Acts
