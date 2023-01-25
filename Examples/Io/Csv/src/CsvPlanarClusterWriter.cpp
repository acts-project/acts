// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvPlanarClusterWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Digitization/DigitizationSourceLink.hpp"
#include "Acts/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvPlanarClusterWriter::CsvPlanarClusterWriter(
    const ActsExamples::CsvPlanarClusterWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputClusters, "CsvPlanarClusterWriter", level),
      m_cfg(config) {
  // inputClusters is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
}

ActsExamples::ProcessCode ActsExamples::CsvPlanarClusterWriter::writeT(
    const AlgorithmContext& ctx,
    const ActsExamples::GeometryIdMultimap<Acts::PlanarModuleCluster>&
        clusters) {
  // retrieve simulated hits
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);

  // open per-event file for all components
  std::string pathHits =
      perEventFilepath(m_cfg.outputDir, "hits.csv", ctx.eventNumber);
  std::string pathCells =
      perEventFilepath(m_cfg.outputDir, "cells.csv", ctx.eventNumber);
  std::string pathTruth =
      perEventFilepath(m_cfg.outputDir, "truth.csv", ctx.eventNumber);

  dfe::NamedTupleCsvWriter<HitData> writerHits(pathHits, m_cfg.outputPrecision);
  dfe::NamedTupleCsvWriter<CellData> writerCells(pathCells,
                                                 m_cfg.outputPrecision);
  dfe::NamedTupleCsvWriter<TruthHitData> writerTruth(pathTruth,
                                                     m_cfg.outputPrecision);

  HitData hit;
  CellData cell;
  TruthHitData truth;
  // will be reused as hit counter
  hit.hit_id = 0;

  for (auto [moduleGeoId, moduleClusters] : groupByModule(clusters)) {
    const Acts::Surface* surfacePtr =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
    if (surfacePtr == nullptr) {
      ACTS_ERROR("Could not find surface for " << moduleGeoId);
      return ProcessCode::ABORT;
    }

    // encoded geometry identifier. same for all hits on the module
    hit.geometry_id = moduleGeoId.value();

    for (const auto& entry : moduleClusters) {
      const Acts::PlanarModuleCluster& cluster = entry.second;
      // local cluster information
      const auto& parameters = cluster.parameters();
      Acts::Vector2 localPos(parameters[0], parameters[1]);
      Acts::Vector3 globalFakeMom(1, 1, 1);
      Acts::Vector3 globalPos =
          surfacePtr->localToGlobal(ctx.geoContext, localPos, globalFakeMom);

      // write global hit information
      hit.x = globalPos.x() / Acts::UnitConstants::mm;
      hit.y = globalPos.y() / Acts::UnitConstants::mm;
      hit.z = globalPos.z() / Acts::UnitConstants::mm;
      hit.t = parameters[2] / Acts::UnitConstants::ns;
      writerHits.append(hit);

      // write local cell information
      cell.hit_id = hit.hit_id;
      for (auto& c : cluster.digitizationCells()) {
        cell.channel0 = c.channel0;
        cell.channel1 = c.channel1;
        // TODO store digitial timestamp once added to the cell definition
        cell.timestamp = 0;
        cell.value = c.data;
        writerCells.append(cell);
      }

      // write hit-particle truth association
      // each hit can have multiple particles, e.g. in a dense environment
      truth.hit_id = hit.hit_id;
      truth.geometry_id = hit.geometry_id;
      const auto& sl = cluster.sourceLink().get<Acts::DigitizationSourceLink>();
      for (auto idx : sl.indices()) {
        auto it = simHits.nth(idx);
        if (it == simHits.end()) {
          ACTS_FATAL("Simulation hit with index " << idx << " does not exist");
          return ProcessCode::ABORT;
        }

        const auto& simHit = *it;
        truth.particle_id = simHit.particleId().value();
        // hit position
        truth.tx = simHit.position().x() / Acts::UnitConstants::mm;
        truth.ty = simHit.position().y() / Acts::UnitConstants::mm;
        truth.tz = simHit.position().z() / Acts::UnitConstants::mm;
        truth.tt = simHit.time() / Acts::UnitConstants::ns;
        // particle four-momentum before interaction
        truth.tpx = simHit.momentum4Before().x() / Acts::UnitConstants::GeV;
        truth.tpy = simHit.momentum4Before().y() / Acts::UnitConstants::GeV;
        truth.tpz = simHit.momentum4Before().z() / Acts::UnitConstants::GeV;
        truth.te = simHit.momentum4Before().w() / Acts::UnitConstants::GeV;
        // particle four-momentum change due to interaction
        const auto delta4 = simHit.momentum4After() - simHit.momentum4Before();
        truth.deltapx = delta4.x() / Acts::UnitConstants::GeV;
        truth.deltapy = delta4.y() / Acts::UnitConstants::GeV;
        truth.deltapz = delta4.z() / Acts::UnitConstants::GeV;
        truth.deltae = delta4.w() / Acts::UnitConstants::GeV;
        // TODO write hit index along the particle trajectory
        truth.index = simHit.index();
        writerTruth.append(truth);
      }

      // increase hit id for next iteration
      hit.hit_id += 1;
    }
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
