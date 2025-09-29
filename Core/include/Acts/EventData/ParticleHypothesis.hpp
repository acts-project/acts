// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/GenericParticleHypothesis.hpp"

namespace Acts {

// TODO In principle the factory methods could provide a reference to a static
// instance which would avoid copying the particle hypothesis and potentially
// save some memory. But constexpr+static seems to require C++2b extension.

/// Specialized particle hypothesis for singly charged particles.
///
/// @note This serves as a factory for common singly charge particles.
class SinglyChargedParticleHypothesis
    : public GenericParticleHypothesis<SinglyCharged> {
 public:
  /// Constructor with explicit mass
  /// @param absPdg The absolute PDG particle code
  /// @param mass The particle mass
  constexpr SinglyChargedParticleHypothesis(PdgParticle absPdg, float mass)
      : GenericParticleHypothesis(absPdg, mass, {}) {}

  /// Constructor with PDG particle code (mass from particle data table)
  /// @param absPdg The absolute PDG particle code
  explicit SinglyChargedParticleHypothesis(PdgParticle absPdg)
      : GenericParticleHypothesis(absPdg) {}

  /// Convert from another particle hypothesis with different charge type
  /// @param other The source particle hypothesis to convert from
  template <typename other_charge_t>
  explicit constexpr SinglyChargedParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  /// Create a muon particle hypothesis
  /// @return Singly charged muon particle hypothesis
  static SinglyChargedParticleHypothesis muon() {
    static const SinglyChargedParticleHypothesis cache(PdgParticle::eMuon);
    return cache;
  }
  /// Create a charged pion particle hypothesis
  /// @return Singly charged pion particle hypothesis
  static SinglyChargedParticleHypothesis pion() {
    static const SinglyChargedParticleHypothesis cache(PdgParticle::ePionPlus);
    return cache;
  }
  /// Create an electron particle hypothesis
  /// @return Singly charged electron particle hypothesis
  static SinglyChargedParticleHypothesis electron() {
    static const SinglyChargedParticleHypothesis cache(PdgParticle::eElectron);
    return cache;
  }
  /// Create a charged kaon particle hypothesis
  /// @return Singly charged kaon particle hypothesis
  static SinglyChargedParticleHypothesis kaon() {
    static const SinglyChargedParticleHypothesis cache(PdgParticle::eKaonPlus);
    return cache;
  }
  /// Create a proton particle hypothesis
  /// @return Singly charged proton particle hypothesis
  static SinglyChargedParticleHypothesis proton() {
    static const SinglyChargedParticleHypothesis cache(PdgParticle::eProton);
    return cache;
  }

  /// Create a charged geantino particle hypothesis (massless charged particle)
  /// @return Singly charged geantino particle hypothesis
  static SinglyChargedParticleHypothesis chargedGeantino() {
    static const SinglyChargedParticleHypothesis cache(PdgParticle::eInvalid,
                                                       0);
    return cache;
  }
};

/// Specialized particle hypothesis for neutral particles.
///
/// @note This serves as a factory for common neutral particles.
class NeutralParticleHypothesis : public GenericParticleHypothesis<Neutral> {
 public:
  /// Constructor with explicit mass
  /// @param absPdg The absolute PDG particle code
  /// @param mass The particle mass
  constexpr NeutralParticleHypothesis(PdgParticle absPdg, float mass)
      : GenericParticleHypothesis(absPdg, mass, {}) {}
  /// Constructor with PDG particle code (mass from particle data table)
  /// @param absPdg The absolute PDG particle code
  explicit NeutralParticleHypothesis(PdgParticle absPdg)
      : GenericParticleHypothesis(absPdg) {}

  /// Convert from another particle hypothesis with different charge type
  /// @param other The source particle hypothesis to convert from
  template <typename other_charge_t>
  explicit constexpr NeutralParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  /// Create a photon particle hypothesis
  /// @return Neutral photon particle hypothesis
  static NeutralParticleHypothesis photon() {
    static const NeutralParticleHypothesis cache(PdgParticle::eGamma);
    return cache;
  }
  /// Create a neutral pion particle hypothesis
  /// @return Neutral pion particle hypothesis
  static NeutralParticleHypothesis pion0() {
    static const NeutralParticleHypothesis cache(PdgParticle::ePionZero);
    return cache;
  }

  /// Create a neutral geantino particle hypothesis (massless neutral particle)
  /// @return Neutral geantino particle hypothesis
  static NeutralParticleHypothesis geantino() {
    static const NeutralParticleHypothesis cache(PdgParticle::eInvalid, 0);
    return cache;
  }
};

/// Specialized particle hypothesis for non-neutral particles.
///
/// @note This serves as a factory for common non-neutral particles.
class NonNeutralChargedParticleHypothesis
    : public GenericParticleHypothesis<NonNeutralCharge> {
 public:
  /// Constructor with explicit mass and charge
  /// @param absPdg The absolute PDG particle code
  /// @param mass The particle mass
  /// @param chargeType The non-neutral charge type
  constexpr NonNeutralChargedParticleHypothesis(PdgParticle absPdg, float mass,
                                                NonNeutralCharge chargeType)
      : GenericParticleHypothesis(absPdg, mass, chargeType) {}
  /// Constructor with PDG particle code (mass from particle data table)
  /// @param absPdg The absolute PDG particle code
  explicit NonNeutralChargedParticleHypothesis(PdgParticle absPdg)
      : GenericParticleHypothesis(absPdg) {}

  /// Convert from another particle hypothesis with different charge type
  /// @param other The source particle hypothesis to convert from
  template <typename other_charge_t>
  explicit constexpr NonNeutralChargedParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  /// Create a muon particle hypothesis
  /// @return Non-neutral charged muon particle hypothesis
  static NonNeutralChargedParticleHypothesis muon() {
    return NonNeutralChargedParticleHypothesis{
        SinglyChargedParticleHypothesis::muon()};
  }
  /// Create a charged pion particle hypothesis
  /// @return Non-neutral charged pion particle hypothesis
  static NonNeutralChargedParticleHypothesis pion() {
    return NonNeutralChargedParticleHypothesis{
        SinglyChargedParticleHypothesis::pion()};
  }
  /// Create an electron particle hypothesis
  /// @return Non-neutral charged electron particle hypothesis
  static NonNeutralChargedParticleHypothesis electron() {
    return NonNeutralChargedParticleHypothesis{
        SinglyChargedParticleHypothesis::electron()};
  }
  /// Create a charged kaon particle hypothesis
  /// @return Non-neutral charged kaon particle hypothesis
  static NonNeutralChargedParticleHypothesis kaon() {
    return NonNeutralChargedParticleHypothesis{
        SinglyChargedParticleHypothesis::kaon()};
  }
  /// Create a proton particle hypothesis
  /// @return Non-neutral charged proton particle hypothesis
  static NonNeutralChargedParticleHypothesis proton() {
    return NonNeutralChargedParticleHypothesis{
        SinglyChargedParticleHypothesis::proton()};
  }

  /// Create a pion-like particle hypothesis with custom charge
  /// @param absQ The absolute charge value
  /// @return Non-neutral charged pion-like particle hypothesis
  static NonNeutralChargedParticleHypothesis pionLike(float absQ) {
    return NonNeutralChargedParticleHypothesis(
        pion().absolutePdg(), pion().mass(), NonNeutralCharge{absQ});
  }

  /// Create a charged geantino particle hypothesis with unit charge
  /// @return Non-neutral charged geantino particle hypothesis
  static NonNeutralChargedParticleHypothesis chargedGeantino() {
    static const auto cache = chargedGeantino(Acts::UnitConstants::e);
    return cache;
  }
  /// Create a charged geantino particle hypothesis with custom charge
  /// @param absQ The absolute charge value
  /// @return Non-neutral charged geantino particle hypothesis
  static NonNeutralChargedParticleHypothesis chargedGeantino(float absQ) {
    return NonNeutralChargedParticleHypothesis(PdgParticle::eInvalid, 0,
                                               NonNeutralCharge{absQ});
  }
};

/// Specialized particle hypothesis for any kind of charged particles.
///
/// @note This serves as a factory for common particles with any kind of charge.
class ParticleHypothesis : public GenericParticleHypothesis<AnyCharge> {
 public:
  /// Constructor with explicit mass and charge
  /// @param absPdg The absolute PDG particle code
  /// @param mass The particle mass
  /// @param chargeType The charge type (any charge)
  constexpr ParticleHypothesis(PdgParticle absPdg, float mass,
                               AnyCharge chargeType)
      : GenericParticleHypothesis(absPdg, mass, chargeType) {}
  /// Constructor with PDG particle code (mass from particle data table)
  /// @param absPdg The absolute PDG particle code
  explicit ParticleHypothesis(PdgParticle absPdg)
      : GenericParticleHypothesis(absPdg) {}

  /// Convert from another particle hypothesis with different charge type
  /// @param other The source particle hypothesis to convert from
  template <typename other_charge_t>
  explicit constexpr ParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  /// Create a muon particle hypothesis
  /// @return Muon particle hypothesis with any charge type
  static ParticleHypothesis muon() {
    return ParticleHypothesis{SinglyChargedParticleHypothesis::muon()};
  }
  /// Create a charged pion particle hypothesis
  /// @return Charged pion particle hypothesis with any charge type
  static ParticleHypothesis pion() {
    return ParticleHypothesis{SinglyChargedParticleHypothesis::pion()};
  }
  /// Create an electron particle hypothesis
  /// @return Electron particle hypothesis with any charge type
  static ParticleHypothesis electron() {
    return ParticleHypothesis{SinglyChargedParticleHypothesis::electron()};
  }
  /// Create a charged kaon particle hypothesis
  /// @return Charged kaon particle hypothesis with any charge type
  static ParticleHypothesis kaon() {
    return ParticleHypothesis{SinglyChargedParticleHypothesis::kaon()};
  }
  /// Create a proton particle hypothesis
  /// @return Proton particle hypothesis with any charge type
  static ParticleHypothesis proton() {
    return ParticleHypothesis{SinglyChargedParticleHypothesis::proton()};
  }

  /// Create a photon particle hypothesis
  /// @return Photon particle hypothesis with any charge type
  static ParticleHypothesis photon() {
    return ParticleHypothesis{NeutralParticleHypothesis::photon()};
  }
  /// Create a neutral pion particle hypothesis
  /// @return Neutral pion particle hypothesis with any charge type
  static ParticleHypothesis pion0() {
    return ParticleHypothesis{NeutralParticleHypothesis::pion0()};
  }

  /// Create a pion-like particle hypothesis with custom charge
  /// @param absQ The absolute charge value
  /// @return Pion-like particle hypothesis with any charge type
  static ParticleHypothesis pionLike(float absQ) {
    return ParticleHypothesis(pion().absolutePdg(), pion().mass(),
                              AnyCharge{absQ});
  }

  /// Create a neutral geantino particle hypothesis (massless neutral particle)
  /// @return Neutral geantino particle hypothesis with any charge type
  static ParticleHypothesis geantino() {
    return ParticleHypothesis{NeutralParticleHypothesis::geantino()};
  }
  /// Create a charged geantino particle hypothesis with unit charge
  /// @return Charged geantino particle hypothesis with any charge type
  static ParticleHypothesis chargedGeantino() {
    static const auto cache = chargedGeantino(Acts::UnitConstants::e);
    return cache;
  }
  /// Create a charged geantino particle hypothesis with custom charge
  /// @param absQ The absolute charge value
  /// @return Charged geantino particle hypothesis with any charge type
  static ParticleHypothesis chargedGeantino(float absQ) {
    return ParticleHypothesis(PdgParticle::eInvalid, 0, AnyCharge{absQ});
  }
};

}  // namespace Acts
