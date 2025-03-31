// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/ParticleData.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <limits>
#include <optional>
#include <ostream>
#include <type_traits>

#include "ParticleDataTable.hpp"

namespace {

static constexpr float extractCharge(float value) {
  // convert three charge to regular charge in native units
  return (value / 3.0f) * Acts::UnitConstants::e;
}

static constexpr float extractMass(float value) {
  return value * Acts::UnitConstants::MeV;
}

}  // namespace

std::optional<float> Acts::findCharge(Acts::PdgParticle pdg) {
  const auto it = kParticlesMapThreeCharge.find(pdg);
  if (it == kParticlesMapThreeCharge.end()) {
    return std::nullopt;
  }
  const auto charge =
      static_cast<std::int32_t>(kParticlesMapThreeCharge.at(pdg));
  return extractCharge(charge);
}

std::optional<float> Acts::findMass(Acts::PdgParticle pdg) {
  const auto it = kParticlesMapMassMeV.find(pdg);
  if (it == kParticlesMapMassMeV.end()) {
    return std::nullopt;
  }
  const auto mass = static_cast<float>(kParticlesMapMassMeV.at(pdg));
  return extractMass(mass);
}

std::optional<std::string_view> Acts::findName(Acts::PdgParticle pdg) {
  const auto it = kParticlesMapName.find(pdg);
  if (it == kParticlesMapName.end()) {
    return std::nullopt;
  }
  return kParticlesMapName.at(pdg);
}

std::optional<Acts::ParticleData> Acts::findParticleData(PdgParticle pdg) {
  const auto it = kParticlesMapThreeCharge.find(pdg);
  if (it == kParticlesMapThreeCharge.end()) {
    return std::nullopt;
  }
  const auto charge = static_cast<float>(kParticlesMapThreeCharge.at(pdg));
  const auto mass = static_cast<float>(kParticlesMapMassMeV.at(pdg));
  const auto name = kParticlesMapName.at(pdg);

  return ParticleData{extractCharge(charge), extractMass(mass), name};
}

std::ostream& Acts::operator<<(std::ostream& os, Acts::PdgParticle pdg) {
  const auto name = Acts::findName(pdg);
  if (name) {
    os << *name;
  } else {
    os << static_cast<std::int32_t>(pdg);
  }
  return os;
}

std::optional<std::string_view> Acts::pdgToShortAbsString(PdgParticle pdg) {
  pdg = makeAbsolutePdgParticle(pdg);
  if (pdg == eElectron) {
    return "e";
  }
  if (pdg == eMuon) {
    return "mu";
  }
  if (pdg == eTau) {
    return "t";
  }
  if (pdg == eGamma) {
    return "g";
  }
  if (pdg == ePionZero) {
    return "pi0";
  }
  if (pdg == ePionPlus) {
    return "pi";
  }
  if (pdg == eKaonPlus) {
    return "K";
  }
  if (pdg == eNeutron) {
    return "n";
  }
  if (pdg == eProton) {
    return "p";
  }
  if (pdg == eLead) {
    return "lead";
  }
  return std::nullopt;
}
