// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/EventGenerator.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <string>

namespace ActsFatras::Synthetic {

namespace {

/// Create the dynamic columns the generator fills, unless they are already
/// there: `SpacePointContainer::clear` keeps them, and `createColumn` throws.
struct Columns {
  Acts::MutableSpacePointColumnProxy<std::uint32_t> layer;
  Acts::MutableSpacePointColumnProxy<float> clusterWidth;
  Acts::MutableSpacePointColumnProxy<float> localPositionY;
  Acts::MutableSpacePointColumnProxy<std::uint32_t> particle;
};

template <typename T>
Acts::MutableSpacePointColumnProxy<T> ensureColumn(
    Acts::SpacePointContainer& spacePoints, const std::string& name) {
  if (spacePoints.hasColumn(name)) {
    return spacePoints.column<T>(name);
  }
  return spacePoints.createColumn<T>(name);
}

Columns ensureColumns(Acts::SpacePointContainer& spacePoints) {
  return Columns{ensureColumn<std::uint32_t>(spacePoints, "layerId"),
                 ensureColumn<float>(spacePoints, "clusterWidth"),
                 ensureColumn<float>(spacePoints, "localPositionY"),
                 ensureColumn<std::uint32_t>(spacePoints, "particleId")};
}

}  // namespace

EventGenerator::EventGenerator(const DetectorLayout& layout,
                               const EventConfig& config)
    : m_layout(&layout),
      m_cfg(config),
      m_primaries(config),
      m_propagator(layout, config),
      m_rng(config.seed) {}

void EventGenerator::reset(const std::uint32_t seed) {
  m_rng.seed(seed);
}

Event EventGenerator::generate() {
  Event event;
  generate(event);
  return event;
}

void EventGenerator::generate(Event& event) {
  constexpr float pi = std::numbers::pi_v<float>;

  const DetectorLayout& layout = *m_layout;
  const MeasurementConfig& cfg = m_cfg.simulation.measurement;

  m_scratch.tracks.clear();
  event.particles.clear();
  m_primaries.generate(m_rng, m_scratch.tracks, event.particles);
  m_propagator.propagate(m_rng, m_scratch, event.particles);

  Acts::SpacePointContainer& spacePoints = event.spacePoints;
  spacePoints.clear();
  spacePoints.createColumns(
      Acts::SpacePointColumns::CopiedFromIndex | Acts::SpacePointColumns::X |
      Acts::SpacePointColumns::Y | Acts::SpacePointColumns::Z |
      Acts::SpacePointColumns::R | Acts::SpacePointColumns::Phi |
      Acts::SpacePointColumns::VarianceZ | Acts::SpacePointColumns::VarianceR);
  // not const: `SpacePointProxy::extra` takes a mutable column proxy by
  // non-const reference
  Columns columns = ensureColumns(spacePoints);
  spacePoints.reserve(static_cast<std::uint32_t>(m_scratch.hits.size()));

  const float measuredVariance = cfg.positionSmearing * cfg.positionSmearing;
  for (const detail::SmearedHit& hit : m_scratch.hits) {
    const bool cylinder =
        layout.layers[hit.layer].shape == SurfaceShape::Cylinder;
    // the helix azimuth is not wrapped while it is propagated
    const float phi = std::remainder(hit.phi, 2.f * pi);
    auto sp = spacePoints.createSpacePoint();
    sp.x() = hit.r * std::cos(phi);
    sp.y() = hit.r * std::sin(phi);
    sp.z() = hit.z;
    sp.r() = hit.r;
    sp.phi() = phi;
    // the hit sits at the nominal position along the normal, so that direction
    // carries no error: a cylinder measures z, a disc r
    sp.varianceZ() = cylinder ? measuredVariance : 0.f;
    sp.varianceR() = cylinder ? 0.f : measuredVariance;
    sp.copiedFromIndex() = sp.index();
    sp.extra(columns.layer) = hit.layer;
    // GBTS reads both; nothing here resolves a cluster, so they stay at zero
    sp.extra(columns.clusterWidth) = 0.f;
    sp.extra(columns.localPositionY) = 0.f;
    sp.extra(columns.particle) = hit.particle;
  }
}

}  // namespace ActsFatras::Synthetic

namespace ActsFatras {

Synthetic::Event Synthetic::generateEvent(const DetectorLayout& layout,
                                          const EventConfig& config) {
  return EventGenerator(layout, config).generate();
}

}  // namespace ActsFatras
