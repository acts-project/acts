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
  /// @return Muon particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis muon() {
    static const ParticleHypothesis cache(PdgParticle::eMuon);
    return cache;
  }
  /// Create a charged pion particle hypothesis
  /// @return Charged pion particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis pion() {
    static const ParticleHypothesis cache(PdgParticle::ePionPlus);
    return cache;
  }
  /// Create an electron particle hypothesis
  /// @return Electron particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis electron() {
    static const ParticleHypothesis cache(PdgParticle::eElectron);
    return cache;
  }
  /// Create a charged kaon particle hypothesis
  /// @return Charged kaon particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis kaon() {
    static const ParticleHypothesis cache(PdgParticle::eKaonPlus);
    return cache;
  }
  /// Create a proton particle hypothesis
  /// @return Proton particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis proton() {
    static const ParticleHypothesis cache(PdgParticle::eProton);
    return cache;
  }

  /// Create a photon particle hypothesis
  /// @return Photon particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis photon() {
    static const ParticleHypothesis cache(PdgParticle::eGamma);
    return cache;
  }
  /// Create a neutral pion particle hypothesis
  /// @return Neutral pion particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis pion0() {
    static const ParticleHypothesis cache(PdgParticle::ePionZero);
    return cache;
  }

  /// Create a pion-like particle hypothesis with custom charge
  /// @param absoluteCharge The absolute charge value
  /// @return Pion-like particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis pionLike(float absoluteCharge) {
    return pion().withAlteredAbsoluteCharge(absoluteCharge);
  }

  /// Create a neutral geantino particle hypothesis (massless neutral particle)
  /// @return Neutral geantino particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis geantino() {
    static const ParticleHypothesis cache(PdgParticle::eInvalid, 0,
                                          ChargeHypothesis{0});
    return cache;
  }
  /// Create a charged geantino particle hypothesis with unit charge
  /// @return Charged geantino particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis chargedGeantino() {
    return geantino().withAlteredAbsoluteCharge(1 * UnitConstants::e);
  }
  /// Create a charged geantino particle hypothesis with custom charge
  /// @param absoluteCharge The absolute charge value
  /// @return Charged geantino particle hypothesis with any charge type
  [[nodiscard]] static ParticleHypothesis chargedGeantino(
      float absoluteCharge) {
    return geantino().withAlteredAbsoluteCharge(absoluteCharge);
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

  /// Create a new particle hypothesis with the same mass and charge but a
  /// different absolute PDG.
  /// @param absPdg The new absolute PDG value
  /// @return A new ParticleHypothesis with the updated absolute PDG
  [[nodiscard]] ParticleHypothesis withAlteredPdg(PdgParticle absPdg) const {
    return ParticleHypothesis(absPdg, mass(), absoluteCharge(), m_momentum);
  }

  /// Create a new particle hypothesis with the same absolute PDG and charge but
  /// a different mass.
  /// @param mass The new mass value
  /// @return A new ParticleHypothesis with the updated mass
  [[nodiscard]] ParticleHypothesis withAlteredMass(float mass) const {
    return ParticleHypothesis(absolutePdg(), mass, absoluteCharge(),
                              m_momentum);
  }

  /// Create a new particle hypothesis with the same absolute PDG and mass but a
  /// different charge.
  /// @param absoluteCharge The new absolute charge value
  /// @return A new ParticleHypothesis with the updated charge
  [[nodiscard]] ParticleHypothesis withAlteredAbsoluteCharge(
      float absoluteCharge) const {
    return ParticleHypothesis(absolutePdg(), mass(), absoluteCharge,
                              m_momentum);
  }

  /// Create a new particle hypothesis with the same absolute PDG and mass but a
  /// different momentum hypothesis.
  /// @param momentum The new momentum hypothesis value
  /// @return A new ParticleHypothesis with the updated momentum
  [[nodiscard]] ParticleHypothesis withMomentumHypothesis(
      double momentum) const {
    return ParticleHypothesis(absolutePdg(), mass(), absoluteCharge(),
                              momentum);
  }

  /// Create a new particle hypothesis with the same absolute PDG and mass but a
  /// different momentum hypothesis.
  /// @param momentum The new optional momentum hypothesis value
  /// @return A new ParticleHypothesis with the updated momentum
  [[nodiscard]] ParticleHypothesis withMomentumHypothesis(
      std::optional<double> momentum) const {
    return ParticleHypothesis(absolutePdg(), mass(), absoluteCharge(),
                              momentum);
  }

  /// Create a new particle hypothesis with the same absolute PDG and mass but
  /// no momentum hypothesis.
  /// @return A new ParticleHypothesis with no momentum hypothesis
  [[nodiscard]] ParticleHypothesis withoutMomentumHypothesis() const {
    return ParticleHypothesis(absolutePdg(), mass(), absoluteCharge(),
                              std::nullopt);
  }

  /// Get the hypothesized absolute PDG.
  /// @return The absolute PDG particle identifier
  [[nodiscard]] constexpr PdgParticle absolutePdg() const noexcept {
    return m_absPdg;
  }

  /// Get the hypothesized mass.
  /// @return The particle mass
  [[nodiscard]] constexpr float mass() const noexcept { return m_mass; }

  /// Get the hypothesized absolute charge.
  /// @return The absolute charge magnitude
  [[nodiscard]] float absoluteCharge() const noexcept {
    return m_charge.absoluteCharge();
  }

  /// Check if the particle hypothesis has a hypothesized momentum.
  /// @return True if a momentum hypothesis is present, false otherwise
  [[nodiscard]] bool hasMomentumHypothesis() const noexcept {
    return m_momentum.has_value();
  }

  /// Get the hypothesized momentum.
  /// @return The particle momentum
  [[nodiscard]] double momentumHypothesis() const { return m_momentum.value(); }

  /// Extracts the signed charge from the `q over p` track parameter using the
  /// charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted signed charge
  [[nodiscard]] constexpr float extractCharge(double qOverP) const noexcept {
    return m_charge.extractCharge(qOverP);
  }

  /// Extracts the particle momentum from the `q over p` track parameter using
  /// the charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted absolute momentum
  [[nodiscard]] constexpr double extractMomentum(double qOverP) const noexcept {
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
  [[nodiscard]] constexpr double qOverP(double momentum,
                                        float signedQ) const noexcept {
    return m_charge.qOverP(momentum, signedQ);
  }

  /// Get the hypothesized charge.
  /// @return Reference to the charge type object
  [[nodiscard]] constexpr const ChargeHypothesis& charge() const noexcept {
    return m_charge;
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
    os << ", mass=" << mass() << ", absCharge=" << absoluteCharge();
    if (hasMomentumHypothesis()) {
      os << ", momentum=" << momentumHypothesis();
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
