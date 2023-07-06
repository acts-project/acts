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

// TODO In principle the factory methods could provide a reference to a static
// instance which would avoid copying the particle hypothesis and potentially
// save some memory. But constexpr+static seems to require C++2b extension.

class SinglyChargedParticleHypothesis
    : public GenericParticleHypothesis<SinglyCharged> {
 public:
  constexpr SinglyChargedParticleHypothesis(PdgParticle absPdg, float mass)
      : GenericParticleHypothesis(absPdg, mass, {}) {}

  template <typename other_charge_t>
  constexpr SinglyChargedParticleHypothesis(
      const GenericParticleHypothesis<other_charge_t>& other)
      : GenericParticleHypothesis(other) {}

  constexpr static SinglyChargedParticleHypothesis muon() {
    return SinglyChargedParticleHypothesis(PdgParticle::eMuon,
                                           105.6583755 * UnitConstants::MeV);
  }
  constexpr static SinglyChargedParticleHypothesis pion() {
    return SinglyChargedParticleHypothesis(PdgParticle::ePionPlus,
                                           139.57039 * UnitConstants::MeV);
  }
  constexpr static SinglyChargedParticleHypothesis electron() {
    return SinglyChargedParticleHypothesis(PdgParticle::ePionPlus,
                                           0.51099895000 * UnitConstants::MeV);
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

  constexpr static NeutralParticleHypothesis photon() {
    return NeutralParticleHypothesis(PdgParticle::eGamma, 0);
  }
  constexpr static NeutralParticleHypothesis pion0() {
    return NeutralParticleHypothesis(PdgParticle::ePionZero,
                                     134.9768 * UnitConstants::MeV);
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

  constexpr static NonNeutralChargedParticleHypothesis muon() {
    return SinglyChargedParticleHypothesis::muon();
  }
  constexpr static NonNeutralChargedParticleHypothesis pion() {
    return SinglyChargedParticleHypothesis::pion();
  }
  constexpr static NonNeutralChargedParticleHypothesis electron() {
    return SinglyChargedParticleHypothesis::electron();
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

  constexpr static ParticleHypothesis muon() {
    return SinglyChargedParticleHypothesis::muon();
  }
  constexpr static ParticleHypothesis pion() {
    return SinglyChargedParticleHypothesis::pion();
  }
  constexpr static ParticleHypothesis electron() {
    return SinglyChargedParticleHypothesis::electron();
  }

  constexpr static ParticleHypothesis photon() {
    return NeutralParticleHypothesis::photon();
  }
  constexpr static ParticleHypothesis pion0() {
    return NeutralParticleHypothesis::pion0();
  }

  constexpr static ParticleHypothesis pionLike(float absQ) {
    return ParticleHypothesis(pion().absPdg(), pion().mass(), absQ);
  }
};

}  // namespace Acts
