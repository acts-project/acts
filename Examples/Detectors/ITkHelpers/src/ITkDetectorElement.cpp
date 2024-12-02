// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/ITkHelpers/ITkDetectorElement.hpp"

#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include <ranges>

namespace ActsExamples {

// Mapping between the barrel-endcap identifier and its unsigned representation
constexpr static std::array<std::pair<unsigned, int>, 3> s_barrelEndcapMap{
    {{0, 0}, {1, 2}, {2, -2}}};

ActsExamples::ITkIdentifier::ITkIdentifier(int hardware, int barrelEndcap,
                                           int layerWheel, int etaModule,
                                           int phiModule, int side) {
  assert((hardware == 0) || (hardware == 1));
  assert((barrelEndcap == 2) || (barrelEndcap == -2) || (barrelEndcap == 0));
  assert(layerWheel >= 0);
  assert(phiModule >= 0);
  assert((side == 0) || (side == 1));

  m_identifier.set(0, hardware);

  auto found = std::ranges::find(s_barrelEndcapMap, barrelEndcap,
                                 &std::pair<unsigned, int>::second);
  if (found == s_barrelEndcapMap.end()) {
    throw std::invalid_argument("Invalid barrel-endcap specifier");
  }
  m_identifier.set(1, found->first);
  m_identifier.set(2, layerWheel);
  m_identifier.set(3, static_cast<std::size_t>(etaModule < 0));
  m_identifier.set(4, std::abs(etaModule));
  m_identifier.set(5, phiModule);
  m_identifier.set(6, side);
}

int ActsExamples::ITkIdentifier::hardware() const {
  return m_identifier.level(0);
}

int ActsExamples::ITkIdentifier::barrelEndcap() const {
  auto found = std::ranges::find(s_barrelEndcapMap, m_identifier.level(1),
                                 &std::pair<unsigned, int>::first);
  if (found == s_barrelEndcapMap.end()) {
    throw std::invalid_argument("Invalid barrel-endcap specifier");
  }
  return found->second;
}

int ActsExamples::ITkIdentifier::layerWheel() const {
  return m_identifier.level(2);
}

int ActsExamples::ITkIdentifier::etaModule() const {
  int sign = (m_identifier.level(3) == 0) ? 1 : -1;
  return sign * m_identifier.level(4);
}

int ActsExamples::ITkIdentifier::phiModule() const {
  return m_identifier.level(5);
}

int ActsExamples::ITkIdentifier::side() const {
  return m_identifier.level(6);
}

std::size_t ActsExamples::ITkIdentifier::value() const {
  return m_identifier.value();
}

std::ostream &operator<<(std::ostream &os, const ITkIdentifier &id) {
  os << "(hw: " << id.hardware() << ", be: " << id.barrelEndcap()
     << ", lw: " << id.layerWheel() << ", em: " << id.etaModule()
     << ", pm: " << id.phiModule() << ", sd: " << id.side() << ")";
  return os;
}

}  // namespace ActsExamples
