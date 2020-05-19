// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "PrintHits.hpp"

#include <vector>

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/IndexContainers.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Range.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Utilities/Logger.hpp"

FW::PrintHits::PrintHits(const FW::PrintHits::Config& cfg,
                         Acts::Logging::Level level)
    : BareAlgorithm("PrintHits", level), m_cfg(cfg) {}

FW::ProcessCode FW::PrintHits::execute(const FW::AlgorithmContext& ctx) const {
  using Clusters = FW::GeometryIdMultimap<Acts::PlanarModuleCluster>;
  using HitParticlesMap = FW::IndexMultimap<ActsFatras::Barcode>;
  using HitIds = std::vector<size_t>;

  const auto& clusters = ctx.eventStore.get<Clusters>(m_cfg.inputClusters);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);
  const auto& hitIds = ctx.eventStore.get<HitIds>(m_cfg.inputHitIds);

  // print hits selected by id
  ACTS_INFO("Hits by id selection")
  size_t hitIdEnd = m_cfg.hitIdStart + m_cfg.hitIdLength;
  for (size_t ihit = m_cfg.hitIdStart; ihit < hitIdEnd; ++ihit) {
    auto hitId = hitIds[ihit];
    auto ic = clusters.nth(ihit);
    if (ic == clusters.end()) {
      break;
    }
    Acts::GeometryID geoId = ic->first;
    const Acts::PlanarModuleCluster& c = ic->second;
    ACTS_INFO("  Cluster " << ihit << " hitId " << hitId << " geoId " << geoId
                           << " size " << c.digitizationCells().size());
    // get all contributing particles
    for (const auto& p : makeRange(hitParticlesMap.equal_range(ihit))) {
      ACTS_INFO("    generating particle " << p.first);
    }
  }

  // print hits within geometry selection
  auto numVolume = selectVolume(clusters, m_cfg.volumeId).size();
  auto numLayer = selectLayer(clusters, m_cfg.volumeId, m_cfg.layerId).size();
  auto rangeModule =
      selectModule(clusters, m_cfg.volumeId, m_cfg.layerId, m_cfg.moduleId);

  ACTS_INFO("Hits total: " << clusters.size());
  ACTS_INFO("Hits in volume " << m_cfg.volumeId << ": " << numVolume);
  ACTS_INFO("Hits in volume " << m_cfg.volumeId << " layer " << m_cfg.layerId
                              << ": " << numLayer);
  ACTS_INFO("Hits in volume " << m_cfg.volumeId << " layer " << m_cfg.layerId
                              << " module " << m_cfg.moduleId << ": "
                              << rangeModule.size());
  // we could also use for (const Acts::PlanarModuleCluster& c : rangeModule)
  // for simplicity, but then we could not get the hit index.
  ACTS_INFO("Hits by geometry selection")
  for (auto ic = rangeModule.begin(); ic != rangeModule.end(); ++ic) {
    auto ihit = clusters.index_of(ic);
    auto hitId = hitIds[ihit];

    Acts::GeometryID geoId = ic->first;
    const Acts::PlanarModuleCluster& c = ic->second;
    ACTS_INFO("  Cluster " << ihit << " hitId " << hitId << " geoId " << geoId
                           << " size " << c.digitizationCells().size());
  }

  return ProcessCode::SUCCESS;
}
