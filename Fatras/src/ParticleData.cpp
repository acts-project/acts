// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Utilities/ParticleData.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <limits>
#include <optional>
#include <ostream>

#include "Acts/Utilities/Units.hpp"

#include "ParticleDataTable.hpp"

// Find an element within a data column using sorted pdg numbers as the index.
template <typename ColumnContainer>
static inline auto findByPdg(int32_t pdg, const ColumnContainer& column)
    -> std::optional<std::decay_t<decltype(column[0])>> {
  // should be a static_assert, but that seems to fail on LLVM
  assert((std::size(column) == kParticlesCount) and "Inconsistent column size");

  auto beg = std::cbegin(kParticlesPdgNumber);
  auto end = std::cend(kParticlesPdgNumber);
  // assumes sorted container of pdg numbers
  auto pos = std::lower_bound(beg, end, pdg);
  if (pos == end) {
    return std::nullopt;
  }
  if (*pos != pdg) {
    return std::nullopt;
  }
  return std::make_optional(column[std::distance(beg, pos)]);
}

float ActsFatras::findCharge(Acts::PdgParticle pdg) {
  const auto q3 = findByPdg(static_cast<int32_t>(pdg), kParticlesThreeCharge);
  if (q3.has_value()) {
    // convert three charge to regular charge in native units
    return (q3.value() / 3.0f) * Acts::UnitConstants::e;
  } else {
    // there is no good default charge. clearly mark the missing value.
    return std::numeric_limits<float>::quiet_NaN();
  }
}

float ActsFatras::findMass(Acts::PdgParticle pdg) {
  const auto mass = findByPdg(static_cast<int32_t>(pdg), kParticlesMassMeV);
  // for medium- to high-pt, zero mass is a reasonable fall-back.
  return mass.value_or(0.0f) * Acts::UnitConstants::MeV;
}

std::string_view ActsFatras::findName(Acts::PdgParticle pdg) {
  return findByPdg(static_cast<int32_t>(pdg), kParticlesName).value_or("");
}

std::ostream& Acts::operator<<(std::ostream& os, Acts::PdgParticle pdg) {
  const auto name = ActsFatras::findName(pdg);
  os << static_cast<int32_t>(pdg);
  if (not name.empty()) {
    os << '|' << name;
  }
  return os;
}
