// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

namespace ActsExamples {

EDM4hepMeasurementWriter::EDM4hepMeasurementWriter(
    const EDM4hepMeasurementWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "EDM4hepMeasurementWriter", level),
      m_cfg(config),
      m_writer(config.outputPath, &m_store) {
  ACTS_VERBOSE("Created output file " << config.outputPath);

  // Input container for measurements is already checked by base constructor
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map input collection");
  }

  m_trackerHitPlaneCollection =
      &m_store.create<edm4hep::TrackerHitPlaneCollection>(
          "ActsTrackerHitsPlane");
  m_writer.registerForWrite("ActsTrackerHitsPlane");

  m_trackerHitRawCollection =
      &m_store.create<edm4hep::TrackerHitCollection>("ActsTrackerHitsRaw");
  m_writer.registerForWrite("ActsTrackerHitsRaw");
}

EDM4hepMeasurementWriter::~EDM4hepMeasurementWriter() {
  m_writer.finish();
}

ActsExamples::ProcessCode EDM4hepMeasurementWriter::endRun() {
  // Write the tree
  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode EDM4hepMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  ClusterContainer clusters;

  if (not m_cfg.inputClusters.empty()) {
    ACTS_VERBOSE("Fetch clusters for writing: " << m_cfg.inputClusters);
    clusters = ctx.eventStore.get<ClusterContainer>(m_cfg.inputClusters);
  }

  ACTS_VERBOSE("Writing " << measurements.size()
                          << " measurements in this event.");

  // TODO we have to distinquish between different events

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& measurement = measurements[hitIdx];

    auto trackerHitPlane = m_trackerHitPlaneCollection->create();

    std::visit(
        [&](const auto& m) {
          Acts::GeometryIdentifier geoId = m.sourceLink().geometryId();

          auto parameters = (m.expander() * m.parameters()).eval();

          trackerHitPlane.setCellID(geoId.value());
          trackerHitPlane.setTime(parameters[Acts::eBoundTime] /
                                  Acts::UnitConstants::ns);
          trackerHitPlane.setU(
              {(float)parameters[Acts::eBoundLoc0], (float)parameters[Acts::eBoundLoc1]});

          // auto covariance = (m.expander() * m.covariance() * m.expander().transpose()).eval();

          if (!clusters.empty()) {
            auto cluster = clusters[hitIdx];

            for (auto& c : cluster.channels) {
              auto trackerHitRaw = m_trackerHitRawCollection->create();
              trackerHitPlane.addToRawHits(trackerHitRaw.getObjectID());

              trackerHitRaw.setCellID(trackerHitPlane.getCellID());
              // don't ask ...
              trackerHitRaw.setType(c.bin[0]);
              trackerHitRaw.setQuality(c.bin[1]);
              trackerHitRaw.setTime(c.activation);
            }
          }
        },
        measurement);
  }

  m_writer.writeEvent();
  m_store.clearCollections();

  return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
