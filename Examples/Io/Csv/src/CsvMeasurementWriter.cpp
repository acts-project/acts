// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "ActsExamples/EventData/AverageSimHits.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <ios>
#include <optional>
#include <stdexcept>

#include <dfe/dfe_io_dsv.hpp>

#include "CsvOutputData.hpp"

ActsExamples::CsvMeasurementWriter::CsvMeasurementWriter(
    const ActsExamples::CsvMeasurementWriter::Config& cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputMeasurements, "CsvMeasurementWriter", lvl), m_cfg(cfg) {
  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map input collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
}

ActsExamples::CsvMeasurementWriter::~CsvMeasurementWriter() {}

ActsExamples::ProcessCode ActsExamples::CsvMeasurementWriter::endRun() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::CsvMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);
  // const auto& hitSimHitsMap = ctx.eventStore.get<IndexMultimap<Index>>(
  //    m_cfg.inputMeasurementSimHitsMap);

  ClusterContainer clusters;
  if (not m_cfg.inputClusters.empty()) {
    clusters = ctx.eventStore.get<ClusterContainer>(m_cfg.inputClusters);
  }

  // Open per-event file for all components
  std::string pathMeasurements =
      perEventFilepath(m_cfg.outputDir, "measurements.csv", ctx.eventNumber);
  std::string pathCells =
      perEventFilepath(m_cfg.outputDir, "cells.csv", ctx.eventNumber);
  std::string pathTruth =
      perEventFilepath(m_cfg.outputDir, "truth.csv", ctx.eventNumber);

  dfe::NamedTupleCsvWriter<MeasurementData> writerMeasurements(
      pathMeasurements, m_cfg.outputPrecision);
  dfe::NamedTupleCsvWriter<CellData> writerCells(pathCells,
                                                 m_cfg.outputPrecision);
  dfe::NamedTupleCsvWriter<TruthHitData> writerTruth(pathTruth,
                                                     m_cfg.outputPrecision);

  MeasurementData meas;
  CellData cell;
  TruthHitData truth;

  // Will be reused as hit counter
  meas.measurement_id = 0;

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& measurement = measurements[hitIdx];

    std::visit(
        [&](const auto& m) {
          Acts::GeometryIdentifier geoId = m.sourceLink().geometryId();
          // Find the corresponding surface
          const Acts::Surface* surfacePtr =
              m_cfg.trackingGeometry->findSurface(geoId);
          if (not surfacePtr) {
            ACTS_ERROR("Could not find surface for " << geoId);
            return;
          }

          // MEASUREMENT information ------------------------------------
          // Encoded geometry identifier. same for all hits on the module
          meas.geometry_id = geoId.value();
          meas.volume_id = geoId.volume();
          meas.layer_id = geoId.layer();
          meas.module_id = geoId.sensitive();

          meas.local_key = 0;
          // Create a full set of parameters
          auto parameters = m.expander() * m.parameters();
          meas.local0 = parameters[Acts::eBoundLoc0];
          meas.local1 = parameters[Acts::eBoundLoc1];
          meas.phi = parameters[Acts::eBoundPhi];
          meas.theta = parameters[Acts::eBoundTheta];
          meas.time = parameters[Acts::eBoundTime] / Acts::UnitConstants::ns;

          auto covariance = m.expander() * m.covariance();
          meas.varLocal0 = covariance(Acts::eBoundLoc0, Acts::eBoundLoc0);
          meas.varLocal1 = covariance(Acts::eBoundLoc1, Acts::eBoundLoc1);
          meas.varPhi = covariance(Acts::eBoundPhi, Acts::eBoundPhi);
          meas.varTheta = covariance(Acts::eBoundTheta, Acts::eBoundTheta);
          meas.varTime = covariance(Acts::eBoundTime, Acts::eBoundTime);
          for (unsigned int ipar = 0;
               ipar < static_cast<unsigned int>(Acts::eBoundSize); ++ipar) {
            if (m.contains(static_cast<Acts::BoundIndices>(ipar))) {
              meas.local_key = ((1 << (ipar + 1)) | meas.local_key);
            }
          }

          writerMeasurements.append(meas);

          // CLUSTER / channel information ------------------------------
          if (not clusters.empty()) {
            auto cluster = clusters[hitIdx];
            cell.hit_id = meas.measurement_id;
            for (auto& c : cluster.channels) {
              cell.channel0 = c.bin[0];
              cell.channel1 = c.bin[1];
              // TODO store digitial timestamp once added to the cell definition
              cell.timestamp = 0;
              cell.value = c.activation;
              writerCells.append(cell);
            }
          }

          // TRUTH information ------------------------------------------
          // @TODO support mulitple associations
          truth.hit_id = meas.measurement_id;
          truth.geometry_id = meas.geometry_id;
          auto idx = m.sourceLink().index();
          auto it = simHits.nth(idx);
          if (it == simHits.end()) {
            ACTS_FATAL("Simulation hit with index " << idx
                                                    << " does not exist");
            return;
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
          const auto delta4 =
              simHit.momentum4After() - simHit.momentum4Before();
          truth.deltapx = delta4.x() / Acts::UnitConstants::GeV;
          truth.deltapy = delta4.y() / Acts::UnitConstants::GeV;
          truth.deltapz = delta4.z() / Acts::UnitConstants::GeV;
          truth.deltae = delta4.w() / Acts::UnitConstants::GeV;
          // @TODO write hit index along the particle trajectory
          truth.index = simHit.index();
          writerTruth.append(truth);

          // Increase counter
          meas.measurement_id += 1;
        },
        measurement);
  }
  return ActsExamples::ProcessCode::SUCCESS;
}
