// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/ParticleData.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <limits>
#include <optional>
#include <ostream>
#include <type_traits>

#include "ParticleDataTable.hpp"

std::optional<float> Acts::findCharge(Acts::PdgParticle pdg) {
  if (auto cached = findCached<float, Type::Charge>(pdg, kParticlesMapCharge);
      cached) {
    return cached;
  }

  const auto it = kParticlesMapCharge.find(pdg);
  if (it == kParticlesMapCharge.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<float> Acts::findMass(Acts::PdgParticle pdg) {
  if (auto cached = findCached<float, Type::Mass>(pdg, kParticlesMapMass);
      cached) {
    return cached;
  }

  const auto it = kParticlesMapMass.find(pdg);
  if (it == kParticlesMapMass.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<std::string_view> Acts::findName(Acts::PdgParticle pdg) {
  if (auto cached =
          findCached<const char* const, Type::Name>(pdg, kParticlesMapName);
      cached) {
    return cached;
  }

  const auto it = kParticlesMapName.find(pdg);
  if (it == kParticlesMapName.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<Acts::ParticleData> Acts::findParticleData(PdgParticle pdg) {
  const auto itCharge = kParticlesMapCharge.find(pdg);
  const auto itMass = kParticlesMapMass.find(pdg);
  const auto itName = kParticlesMapName.find(pdg);

  if (itCharge == kParticlesMapCharge.end() ||
      itMass == kParticlesMapMass.end() || itName == kParticlesMapName.end()) {
    return std::nullopt;
  }
  const auto charge = static_cast<float>(itCharge->second);
  const auto mass = static_cast<float>(itMass->second);
  const auto name = itName->second;

  return ParticleData{charge, mass, name};
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
