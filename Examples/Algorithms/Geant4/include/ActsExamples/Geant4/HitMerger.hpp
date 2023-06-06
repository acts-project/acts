// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include "G4VPhysicalVolume.hh"

#include <numeric>

namespace ActsExamples {

class HitMerger {
  ActsFatras::Barcode m_barcode;
  Acts::GeometryIdentifier m_geoid;
  const G4VPhysicalVolume* m_volume;
  std::vector<ActsFatras::Hit> m_hits;

 public:
  const G4VPhysicalVolume *geant4Volume() const {
      return m_volume;
  }

  const std::vector<ActsFatras::Hit> &hits() const {
      return m_hits;
  }

  ActsFatras::Barcode particleID() const {
    return m_barcode;
  }

  HitMerger(ActsFatras::Barcode barcode, Acts::GeometryIdentifier geoid, const G4VPhysicalVolume *volume)
      : m_barcode(barcode), m_geoid(geoid), m_volume(volume) {}

  void addHit(ActsFatras::Hit &&hit) {
    if (hit.particleId() != m_barcode or hit.geometryId() != m_geoid) {
      throw std::runtime_error("invalid use of hit merger");
    }
    m_hits.emplace_back(hit);
  }

  ActsFatras::Hit mergeHit(int32_t index) const {
    assert(not m_hits.empty());

    auto pos4 =
        std::accumulate(m_hits.begin(), m_hits.end(), Acts::Vector4::Zero().eval(),
                        [](const Acts::Vector4 &p, const ActsFatras::Hit &h) {
                          return p + h.fourPosition();
                        }) /
        m_hits.size();
    return ActsFatras::Hit(m_geoid, m_barcode, pos4,
                           m_hits.front().momentum4Before(),
                           m_hits.back().momentum4After(), index);
  }
};

}  // namespace ActsExamples
