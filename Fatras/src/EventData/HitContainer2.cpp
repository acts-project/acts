// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/HitContainer2.hpp"

namespace ActsFatras {

void HitContainer2::reserve(std::uint32_t size) noexcept {
  m_hits.reserve(size);
}

void HitContainer2::clear() noexcept {
  m_hits.clear();
  m_hitsByParticles.clear();
  m_hitsBySurfaces.clear();
}

Hit& HitContainer2::push_back(const Hit& hit) {
  m_hits.push_back(hit);
  const auto hitIndex = static_cast<HitIndex>(m_hits.size() - 1);
  m_hitsByParticles[hit.particleId2()].push_back(hitIndex);
  m_hitsBySurfaces[hit.geometryId()].push_back(hitIndex);
  return m_hits.back();
}

Hit& HitContainer2::emplace_back(const Acts::GeometryIdentifier geometryId,
                                 Barcode particleBarcode,
                                 const ParticleIndex2 particleId,
                                 const Acts::Vector4& pos4,
                                 const Acts::Vector4& before4,
                                 const Acts::Vector4& after4,
                                 const std::int32_t index) {
  m_hits.emplace_back(geometryId, particleBarcode, particleId, pos4, before4,
                      after4, index);
  const auto hitIndex = static_cast<HitIndex>(m_hits.size() - 1);
  m_hitsByParticles[particleId].push_back(hitIndex);
  m_hitsBySurfaces[geometryId].push_back(hitIndex);
  return m_hits.back();
}

std::span<const HitIndex> HitContainer2::hitIndicesByParticle(
    ParticleIndex2 particleId) const {
  auto it = m_hitsByParticles.find(particleId);
  if (it != m_hitsByParticles.end()) {
    return it->second;
  }
  return {};
}

std::span<const HitIndex> HitContainer2::hitIndicesBySurface(
    Acts::GeometryIdentifier geometryId) const {
  auto it = m_hitsBySurfaces.find(geometryId);
  if (it != m_hitsBySurfaces.end()) {
    return it->second;
  }
  return {};
}

}  // namespace ActsFatras
