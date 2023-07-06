// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/GenericParticleHypothesis.hpp"

namespace Acts {

class SinglyChargedParticleHypothesis
    : public GenericParticleHypothesis<SinglyCharged> {
 public:
  constexpr SinglyChargedParticleHypothesis(PdgParticle absPdg, float mass)
      : GenericParticleHypothesis(absPdg, mass, {}) {}

  template <typename other_charge_t>
  constexpr SinglyChargedParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  constexpr static const SinglyChargedParticleHypothesis& muon() {
    static constexpr auto p = SinglyChargedParticleHypothesis(
        PdgParticle::eMuon, 105.6583755 * UnitConstants::MeV);
    return p;
  }
  constexpr static SinglyChargedParticleHypothesis pion() {
    static constexpr auto p = SinglyChargedParticleHypothesis(
        PdgParticle::ePionPlus, 139.57039 * UnitConstants::MeV);
    return p;
  }
  constexpr static SinglyChargedParticleHypothesis electron() {
    static constexpr auto p = SinglyChargedParticleHypothesis(
        PdgParticle::ePionPlus, 0.51099895000 * UnitConstants::MeV);
    return p;
  }
};

class NeutralParticleHypothesis : public GenericParticleHypothesis<Neutral> {
 public:
  constexpr NeutralParticleHypothesis(PdgParticle absPdg, float mass)
      : GenericParticleHypothesis(absPdg, mass, {}) {}

  template <typename other_charge_t>
  constexpr NeutralParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  constexpr static const NeutralParticleHypothesis& photon() {
    static constexpr auto p = NeutralParticleHypothesis(PdgParticle::eGamma, 0);
    return p;
  }
  constexpr static const NeutralParticleHypothesis& pion0() {
    static constexpr auto p = NeutralParticleHypothesis(
        PdgParticle::ePionZero, 134.9768 * UnitConstants::MeV);
    return p;
  }
};

class NonNeutralChargedParticleHypothesis
    : public GenericParticleHypothesis<NonNeutralCharge> {
 public:
  constexpr NonNeutralChargedParticleHypothesis(PdgParticle absPdg, float mass,
                                                NonNeutralCharge chargeType)
      : GenericParticleHypothesis(absPdg, mass, chargeType) {}

  template <typename other_charge_t>
  constexpr NonNeutralChargedParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  constexpr static const NonNeutralChargedParticleHypothesis& muon() {
    static constexpr NonNeutralChargedParticleHypothesis p =
        SinglyChargedParticleHypothesis::muon();
    return p;
  }
  constexpr static const NonNeutralChargedParticleHypothesis& pion() {
    static constexpr NonNeutralChargedParticleHypothesis p =
        SinglyChargedParticleHypothesis::pion();
    return p;
  }
  constexpr static const NonNeutralChargedParticleHypothesis& electron() {
    static constexpr NonNeutralChargedParticleHypothesis p =
        SinglyChargedParticleHypothesis::electron();
    return p;
  }

  template <int absQ>
  constexpr static const NonNeutralChargedParticleHypothesis& pionLike() {
    static_assert(absQ > 0, "absQ must be bigger than 0 for NonNeutralCharge");
    static constexpr auto p = NonNeutralChargedParticleHypothesis(
        pion().absPdg(), pion().mass(), absQ);
    return p;
  }
  constexpr static NonNeutralChargedParticleHypothesis pionLike(float absQ) {
    return NonNeutralChargedParticleHypothesis(pion().absPdg(), pion().mass(),
                                               absQ);
  }
};

class ParticleHypothesis : public GenericParticleHypothesis<AnyCharge> {
 public:
  constexpr ParticleHypothesis(PdgParticle absPdg, float mass,
                               AnyCharge chargeType)
      : GenericParticleHypothesis(absPdg, mass, chargeType) {}

  template <typename other_charge_t>
  constexpr ParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  constexpr static const ParticleHypothesis& muon() {
    static constexpr ParticleHypothesis p =
        SinglyChargedParticleHypothesis::muon();
    return p;
  }
  constexpr static const ParticleHypothesis& pion() {
    static constexpr ParticleHypothesis p =
        SinglyChargedParticleHypothesis::pion();
    return p;
  }
  constexpr static const ParticleHypothesis& electron() {
    static constexpr ParticleHypothesis p =
        SinglyChargedParticleHypothesis::electron();
    return p;
  }

  constexpr static const ParticleHypothesis& photon() {
    static constexpr ParticleHypothesis p = NeutralParticleHypothesis::photon();
    return p;
  }
  constexpr static const ParticleHypothesis& pion0() {
    static constexpr ParticleHypothesis p = NeutralParticleHypothesis::pion0();
    return p;
  }

  template <int absQ>
  constexpr static const ParticleHypothesis& pionLike() {
    static_assert(absQ >= 0, "absQ must be positive for AnyCharge");
    static constexpr auto p =
        ParticleHypothesis(pion().absPdg(), pion().mass(), absQ);
    return p;
  }
  constexpr static ParticleHypothesis pionLike(float absQ) {
    return ParticleHypothesis(pion().absPdg(), pion().mass(), absQ);
  }
};

}  // namespace Acts
