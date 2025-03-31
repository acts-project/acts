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

namespace {

enum class Type {
  Charge,
  Mass,
  Name,
};

// Template function depending on the PDG and *type* of data
// so we can use statically cached values
template <typename T, Acts::PdgParticle pdg, Type type>
std::optional<T> findCachedImpl(const std::map<std::int32_t, T>& map) {
  const static std::optional<T> value = [&map]() -> std::optional<T> {
    const auto it = map.find(pdg);
    if (it == map.end()) {
      return std::nullopt;
    }
    return it->second;
  }();

  return value;
}

// Cache lookup for particle data
// Uses a switch statement to map the PDG code to the correct cached value
template <typename T, Type type>
std::optional<T> findCached(Acts::PdgParticle pdg,
                            const std::map<std::int32_t, T>& map) {
  using enum Acts::PdgParticle;
  switch (pdg) {
    case eElectron:
      return findCachedImpl<T, eElectron, type>(map);
    case ePositron:
      return findCachedImpl<T, ePositron, type>(map);
    case eMuon:
      return findCachedImpl<T, eMuon, type>(map);
    case eAntiMuon:
      return findCachedImpl<T, eAntiMuon, type>(map);
    case eTau:
      return findCachedImpl<T, eTau, type>(map);
    case eAntiTau:
      return findCachedImpl<T, eAntiTau, type>(map);
    case eGamma:
      return findCachedImpl<T, eGamma, type>(map);
    case ePionZero:
      return findCachedImpl<T, ePionZero, type>(map);
    case ePionPlus:
      return findCachedImpl<T, ePionPlus, type>(map);
    case ePionMinus:
      return findCachedImpl<T, ePionMinus, type>(map);
    case eKaonPlus:
      return findCachedImpl<T, eKaonPlus, type>(map);
    case eKaonMinus:
      return findCachedImpl<T, eKaonMinus, type>(map);
    case eNeutron:
      return findCachedImpl<T, eNeutron, type>(map);
    case eAntiNeutron:
      return findCachedImpl<T, eAntiNeutron, type>(map);
    case eProton:
      return findCachedImpl<T, eProton, type>(map);
    case eAntiProton:
      return findCachedImpl<T, eAntiProton, type>(map);
    case eLead:
      return findCachedImpl<T, eLead, type>(map);
    default:
      return std::nullopt;
  }
}

}  // namespace

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

  return ParticleData{itCharge->second, itMass->second, itName->second};
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
