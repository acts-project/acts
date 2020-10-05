// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "HitsPrinter.hpp"

#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexContainers.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <vector>

ActsExamples::HitsPrinter::HitsPrinter(
    const ActsExamples::HitsPrinter::Config& cfg, Acts::Logging::Level level)
    : BareAlgorithm("HitsPrinter", level), m_cfg(cfg) {
  if (m_cfg.inputClusters.empty()) {
    throw std::invalid_argument("Input clusters collection is not configured");
  }
  if (m_cfg.inputHitParticlesMap.empty()) {
    throw std::invalid_argument(
        "Input hit-particles map collection is not configured");
  }
  if (m_cfg.inputHitIds.empty()) {
    throw std::invalid_argument("Input hit ids collection is not configured");
  }
}

ActsExamples::ProcessCode ActsExamples::HitsPrinter::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using Clusters = ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>;
  using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;
  using HitIds = std::vector<size_t>;

  const auto& clusters = ctx.eventStore.get<Clusters>(m_cfg.inputClusters);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);
  const auto& hitIds = ctx.eventStore.get<HitIds>(m_cfg.inputHitIds);

  if (clusters.size() != hitIds.size()) {
    ACTS_ERROR(
        "event "
        << ctx.eventNumber
        << " input clusters and hit ids collections have inconsistent size");
    return ProcessCode::ABORT;
  }
  ACTS_INFO("event " << ctx.eventNumber << " collection '"
                     << m_cfg.inputClusters << "' contains " << clusters.size()
                     << " hits");

  // print hits selected by index
  if (0 < m_cfg.selectIndexLength) {
    size_t ihit = m_cfg.selectIndexStart;
    size_t nend = std::min(clusters.size(), ihit + m_cfg.selectIndexLength);

    if (nend <= ihit) {
      ACTS_WARNING("event "
                   << ctx.eventNumber << " collection '" << m_cfg.inputClusters
                   << " hit index selection is outside the available range");
    } else {
      ACTS_INFO("event " << ctx.eventNumber << " hits by index selection");

      Acts::GeometryIdentifier prevGeoId;
      for (; ihit < nend; ++ihit) {
        auto hitId = hitIds.at(ihit);
        auto ic = clusters.nth(ihit);
        if (ic == clusters.end()) {
          break;
        }
        Acts::GeometryIdentifier geoId = ic->first;
        const Acts::PlanarModuleCluster& c = ic->second;
        if (geoId != prevGeoId) {
          ACTS_INFO("on geometry id " << geoId);
          prevGeoId = geoId;
        }
        ACTS_INFO("  cluster " << ihit << " hit " << hitId << " size "
                               << c.digitizationCells().size());
        // get all contributing particles
        for (auto [hid, pid] : makeRange(hitParticlesMap.equal_range(ihit))) {
          ACTS_INFO("    generated by particle " << pid);
        }
      }
    }
  }

  // print hits within geometry selection
  auto geoSelection = makeRange(clusters.begin(), clusters.begin());
  if (m_cfg.selectModule) {
    geoSelection = selectModule(clusters, m_cfg.selectVolume, m_cfg.selectLayer,
                                m_cfg.selectModule);
  } else if (m_cfg.selectLayer) {
    geoSelection = selectLayer(clusters, m_cfg.selectVolume, m_cfg.selectLayer);
  } else if (m_cfg.selectVolume) {
    geoSelection = selectVolume(clusters, m_cfg.selectVolume);
  }
  if (not geoSelection.empty()) {
    ACTS_INFO("event " << ctx.eventNumber << " collection '"
                       << m_cfg.inputClusters << "' contains "
                       << geoSelection.size() << " hits in volume "
                       << m_cfg.selectVolume << " layer " << m_cfg.selectLayer
                       << " module " << m_cfg.selectModule);
    // we could also use for (const Acts::PlanarModuleCluster& c : rangeModule)
    // for simplicity, but then we could not get the hit index.
    Acts::GeometryIdentifier prevGeoId;
    for (auto ic = geoSelection.begin(); ic != geoSelection.end(); ++ic) {
      auto ihit = clusters.index_of(ic);
      auto hitId = hitIds[ihit];

      Acts::GeometryIdentifier geoId = ic->first;
      const Acts::PlanarModuleCluster& c = ic->second;
      if (geoId != prevGeoId) {
        ACTS_INFO("on geometry id " << geoId);
        prevGeoId = geoId;
      }
      ACTS_INFO("  cluster " << ihit << " hit " << hitId << " size "
                             << c.digitizationCells().size());
      // get all contributing particles
      for (auto [hid, pid] : makeRange(hitParticlesMap.equal_range(ihit))) {
        ACTS_INFO("    generated by particle " << pid);
      }
    }
  }

  return ProcessCode::SUCCESS;
}
