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
    const ActsExamples::CsvMeasurementWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "CsvMeasurementWriter", level),
      m_cfg(config) {
  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map input collection");
  }

  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
}

ActsExamples::CsvMeasurementWriter::~CsvMeasurementWriter() = default;

ActsExamples::ProcessCode ActsExamples::CsvMeasurementWriter::finalize() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::CsvMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  ClusterContainer clusters;

  // Open per-event file for all components
  std::string pathMeasurements =
      perEventFilepath(m_cfg.outputDir, "measurements.csv", ctx.eventNumber);
  std::string pathMeasurementSimHitMap = perEventFilepath(
      m_cfg.outputDir, "measurement-simhit-map.csv", ctx.eventNumber);

  dfe::NamedTupleCsvWriter<MeasurementData> writerMeasurements(
      pathMeasurements, m_cfg.outputPrecision);

  std::optional<dfe::NamedTupleCsvWriter<CellData>> writerCells{std::nullopt};
  if (not m_cfg.inputClusters.empty()) {
    ACTS_VERBOSE(
        "Set up writing of clusters from collection: " << m_cfg.inputClusters);
    clusters = m_inputClusters(ctx);
    std::string pathCells =
        perEventFilepath(m_cfg.outputDir, "cells.csv", ctx.eventNumber);
    writerCells =
        dfe::NamedTupleCsvWriter<CellData>{pathCells, m_cfg.outputPrecision};
  }

  dfe::NamedTupleCsvWriter<MeasurementSimHitLink> writerMeasurementSimHitMap(
      pathMeasurementSimHitMap, m_cfg.outputPrecision);

  MeasurementData meas;
  CellData cell;

  // Will be reused as hit counter
  meas.measurement_id = 0;

  ACTS_VERBOSE("Writing " << measurements.size()
                          << " measurements in this event.");

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& measurement = measurements[hitIdx];

    auto simHitIndices = makeRange(measurementSimHitsMap.equal_range(hitIdx));
    for (auto [_, simHitIdx] : simHitIndices) {
      writerMeasurementSimHitMap.append({hitIdx, simHitIdx});
    }

    std::visit(
        [&](const auto& m) {
          Acts::GeometryIdentifier geoId = m.sourceLink().geometryId();
          // MEASUREMENT information ------------------------------------

          // Encoded geometry identifier. same for all hits on the module
          meas.geometry_id = geoId.value();
          meas.local_key = 0;
          // Create a full set of parameters
          auto parameters = (m.expander() * m.parameters()).eval();
          meas.local0 = parameters[Acts::eBoundLoc0];
          meas.local1 = parameters[Acts::eBoundLoc1];
          meas.phi = parameters[Acts::eBoundPhi];
          meas.theta = parameters[Acts::eBoundTheta];
          meas.time = parameters[Acts::eBoundTime] / Acts::UnitConstants::ns;

          auto covariance =
              (m.expander() * m.covariance() * m.expander().transpose()).eval();
          meas.var_local0 = covariance(Acts::eBoundLoc0, Acts::eBoundLoc0);
          meas.var_local1 = covariance(Acts::eBoundLoc1, Acts::eBoundLoc1);
          meas.var_phi = covariance(Acts::eBoundPhi, Acts::eBoundPhi);
          meas.var_theta = covariance(Acts::eBoundTheta, Acts::eBoundTheta);
          meas.var_time = covariance(Acts::eBoundTime, Acts::eBoundTime);
          for (unsigned int ipar = 0;
               ipar < static_cast<unsigned int>(Acts::eBoundSize); ++ipar) {
            if (m.contains(static_cast<Acts::BoundIndices>(ipar))) {
              meas.local_key = ((1 << (ipar + 1)) | meas.local_key);
            }
          }

          writerMeasurements.append(meas);

          // CLUSTER / channel information ------------------------------
          if (not clusters.empty() && writerCells) {
            auto cluster = clusters[hitIdx];
            cell.geometry_id = meas.geometry_id;
            cell.hit_id = meas.measurement_id;
            for (auto& c : cluster.channels) {
              cell.channel0 = c.bin[0];
              cell.channel1 = c.bin[1];
              // TODO store digitial timestamp once added to the cell definition
              cell.timestamp = 0;
              cell.value = c.activation;
              writerCells->append(cell);
            }
          }
          // Increase counter
          meas.measurement_id += 1;
        },
        measurement);
  }
  return ActsExamples::ProcessCode::SUCCESS;
}
