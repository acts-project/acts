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

#include <cassert>
#include <cmath>
#include <iostream>

namespace Acts {

/// @ingroup eventdata
/// @defgroup eventdata-charge Charge interpretation for track parameters
///
/// Track parameters store a single coefficient that describes charge and
/// momentum. This is either charge/momentum or 1/momentum, but the
/// interpretation depends on what type of particle is described. In this code
/// base this coefficient is always referred to as `qOverP` (or
/// charge-over-momentum) even for uncharged particles. The following types are
/// used to restrict the particle charge magnitude (at compile time) and support
/// the umambigous extraction of charge and absolute momentum from said track
/// parameter coefficient.
///
/// All types are designed to be interchangeable. Each one can be
/// constructed with the input charge magnitude
///
/// ```cpp
/// Charge c(1_e);
/// ```
///
/// and can then be used to extract the charge value
///
/// ```cpp
/// auto q = c.extractCharge(qOverP);
/// ```
///
/// or the absolute momentum
///
/// ```cpp
/// auto p = c.extractMomentum(qOverP);
/// ```
///
/// from the charge-over-momentum track parameter.
///
/// @{

/// Charge and momentum interpretation for arbitrarily charged particles.
///
/// Only a charge magnitude identical to zero is interpreted as representing a
/// neutral particle. This avoids ambiguities that might arise from using an
/// approximate comparison with an arbitrary epsilon.
class ChargeHypothesis {
 public:
  /// Construct with the magnitude of the input charge.
  /// @param absoluteCharge The absolute value of the charge magnitude
  constexpr explicit ChargeHypothesis(float absoluteCharge) noexcept
      : m_absoluteCharge{absoluteCharge} {
    assert(absoluteCharge >= 0 &&
           "Input charge magnitude must be zero or positive");
  }

  /// Get the absolute charge magnitude
  /// @return Absolute charge magnitude (0 for neutral particles)
  constexpr float absoluteCharge() const noexcept { return m_absoluteCharge; }

  /// Extract the signed charge from q/p
  /// @param qOverP Charge over momentum
  /// @return Signed charge with correct magnitude (0 for neutral)
  constexpr float extractCharge(double qOverP) const noexcept {
    return static_cast<float>(std::copysign(m_absoluteCharge, qOverP));
  }
  /// Extract momentum magnitude from q/p
  /// @param qOverP Charge over momentum
  /// @return Momentum magnitude (handles both charged and neutral particles)
  constexpr double extractMomentum(double qOverP) const noexcept {
    return (m_absoluteCharge != 0.0f) ? extractCharge(qOverP) / qOverP
                                      : 1.0 / qOverP;
  }

  /// Compute q/p from momentum and signed charge
  /// @param momentum Particle momentum magnitude
  /// @param signedQ Signed charge (must match stored charge magnitude)
  /// @return Charge over momentum (handles both charged and neutral particles)
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    assert(std::abs(signedQ) == m_absoluteCharge && "inconsistent charge");
    return (m_absoluteCharge != 0.0f) ? signedQ / momentum : 1.0 / momentum;
  }

  /// Compare for equality.
  /// @param rhs The ChargeHypothesis to compare to
  /// @return True if the two ChargeHypothesis objects are equal, false otherwise
  constexpr bool operator==(const ChargeHypothesis& rhs) const noexcept =
      default;

 private:
  float m_absoluteCharge{};
};

// TODO In principle the factory methods could provide a reference to a static
// instance which would avoid copying the particle hypothesis and potentially
// save some memory. But constexpr+static seems to require C++2b extension.

/// @brief Particle hypothesis used in reconstruction
///
/// The reconstruction hypothesis consists of absolute PDG code, mass and
/// absolute charge.
class ParticleHypothesis {
 public:
  /// Create a muon particle hypothesis
  /// @return Muon particle hypothesis with any charge type
  static ParticleHypothesis muon() {
    static const ParticleHypothesis cache(PdgParticle::eMuon);
    return cache;
  }
  /// Create a charged pion particle hypothesis
  /// @return Charged pion particle hypothesis with any charge type
  static ParticleHypothesis pion() {
    static const ParticleHypothesis cache(PdgParticle::ePionPlus);
    return cache;
  }
  /// Create an electron particle hypothesis
  /// @return Electron particle hypothesis with any charge type
  static ParticleHypothesis electron() {
    static const ParticleHypothesis cache(PdgParticle::eElectron);
    return cache;
  }
  /// Create a charged kaon particle hypothesis
  /// @return Charged kaon particle hypothesis with any charge type
  static ParticleHypothesis kaon() {
    static const ParticleHypothesis cache(PdgParticle::eKaonPlus);
    return cache;
  }
  /// Create a proton particle hypothesis
  /// @return Proton particle hypothesis with any charge type
  static ParticleHypothesis proton() {
    static const ParticleHypothesis cache(PdgParticle::eProton);
    return cache;
  }

  /// Create a photon particle hypothesis
  /// @return Photon particle hypothesis with any charge type
  static ParticleHypothesis photon() {
    static const ParticleHypothesis cache(PdgParticle::eGamma);
    return cache;
  }
  /// Create a neutral pion particle hypothesis
  /// @return Neutral pion particle hypothesis with any charge type
  static ParticleHypothesis pion0() {
    static const ParticleHypothesis cache(PdgParticle::ePionZero);
    return cache;
  }

  /// Create a pion-like particle hypothesis with custom charge
  /// @param absoluteCharge The absolute charge value
  /// @return Pion-like particle hypothesis with any charge type
  static ParticleHypothesis pionLike(float absoluteCharge) {
    return ParticleHypothesis(pion().absolutePdg(), pion().mass(),
                              ChargeHypothesis{absoluteCharge});
  }

  /// Create a neutral geantino particle hypothesis (massless neutral particle)
  /// @return Neutral geantino particle hypothesis with any charge type
  static ParticleHypothesis geantino() {
    ParticleHypothesis cache(PdgParticle::eInvalid, 0, ChargeHypothesis{0});
    return cache;
  }
  /// Create a charged geantino particle hypothesis with unit charge
  /// @return Charged geantino particle hypothesis with any charge type
  static ParticleHypothesis chargedGeantino() {
    static const auto cache = chargedGeantino(Acts::UnitConstants::e);
    return cache;
  }
  /// Create a charged geantino particle hypothesis with custom charge
  /// @param absoluteCharge The absolute charge value
  /// @return Charged geantino particle hypothesis with any charge type
  static ParticleHypothesis chargedGeantino(float absoluteCharge) {
    return ParticleHypothesis(PdgParticle::eInvalid, 0,
                              ChargeHypothesis{absoluteCharge});
  }

  /// Creates a particle hypothesis using absolute PDG, mass and the charge
  /// type.
  ///
  /// @param absPdg the absolute PDG
  /// @param mass the particle mass
  /// @param absCharge the absolute charge
  constexpr ParticleHypothesis(PdgParticle absPdg, float mass, float absCharge)
      : m_absPdg{absPdg}, m_mass{mass}, m_charge{absCharge} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Creates a particle hypothesis using absolute PDG, mass and the charge
  /// type.
  ///
  /// @param absPdg the absolute PDG
  /// @param mass the particle mass
  /// @param charge the charge type
  constexpr ParticleHypothesis(PdgParticle absPdg, float mass,
                               ChargeHypothesis charge)
      : m_absPdg{absPdg}, m_mass{mass}, m_charge{charge} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Creates a particle hypothesis using the absolute PDG.
  /// The mass and charge is looked up using @ref findMass and @ref findCharge.
  /// If the lookup fails an exception is thrown.
  ///
  /// @param absPdg the absolute PDG
  explicit ParticleHypothesis(PdgParticle absPdg)
      : m_absPdg{absPdg},
        m_mass{findMass(absPdg).value()},
        m_charge{std::abs(findCharge(absPdg).value())} {
    assert(absPdg == makeAbsolutePdgParticle(absPdg) &&
           "pdg is expected to be absolute");
  }

  /// Get the hypothesized absolute PDG.
  /// @return The absolute PDG particle identifier
  constexpr PdgParticle absolutePdg() const noexcept { return m_absPdg; }

  /// Get the hypothesized mass.
  /// @return The particle mass in natural units
  constexpr float mass() const noexcept { return m_mass; }

  /// Get the hypothesized absolute charge.
  /// @return The absolute charge magnitude
  float absoluteCharge() const noexcept { return m_charge.absoluteCharge(); }

  /// Extracts the signed charge from the `q over p` track parameter using the
  /// charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted signed charge
  template <typename T>
  constexpr auto extractCharge(T qOverP) const noexcept {
    return m_charge.extractCharge(qOverP);
  }

  /// Extracts the particle momentum from the `q over p` track parameter using
  /// the charge hypothesis.
  ///
  /// @param qOverP the `q over p` track parameter.
  /// @return The extracted absolute momentum
  template <typename T>
  constexpr auto extractMomentum(T qOverP) const noexcept {
    return m_charge.extractMomentum(qOverP);
  }

  /// Calculate the `q over p` track parameter with the given absolute momentum
  /// and charge.
  ///
  /// @param momentum the absolute momentum.
  /// @param signedQ the signed charge.
  /// @return The calculated charge over momentum ratio
  template <typename P, typename Q>
  constexpr auto qOverP(P momentum, Q signedQ) const noexcept {
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
    os << ", mass=" << mass() << ", absCharge=" << absoluteCharge() << "}";
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

  friend bool operator==(const ParticleHypothesis& lhs,
                         const ParticleHypothesis& rhs) {
    return (lhs.m_absPdg == rhs.m_absPdg) && (lhs.m_mass == rhs.m_mass) &&
           (lhs.m_charge == rhs.m_charge);
  }
};

/// @}

}  // namespace Acts
