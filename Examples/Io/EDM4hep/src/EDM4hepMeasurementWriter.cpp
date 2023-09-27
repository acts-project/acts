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
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

namespace ActsExamples {

EDM4hepMeasurementWriter::EDM4hepMeasurementWriter(
    const EDM4hepMeasurementWriter::Config& config, Acts::Logging::Level level)
    : WriterT(config.inputMeasurements, "EDM4hepMeasurementWriter", level),
      m_cfg(config),
      m_writer(config.outputPath, &m_store) {
  ACTS_VERBOSE("Created output file " << config.outputPath);

  // Input container for measurements is already checked by base constructor
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);

  m_trackerHitPlaneCollection =
      &m_store.create<edm4hep::TrackerHitPlaneCollection>(
          "ActsTrackerHitsPlane");
  m_writer.registerForWrite("ActsTrackerHitsPlane");

  m_trackerHitRawCollection =
      &m_store.create<edm4hep::TrackerHitCollection>("ActsTrackerHitsRaw");
  m_writer.registerForWrite("ActsTrackerHitsRaw");
}

ActsExamples::ProcessCode EDM4hepMeasurementWriter::finalize() {
  m_writer.finish();

  return ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode EDM4hepMeasurementWriter::writeT(
    const AlgorithmContext& ctx, const MeasurementContainer& measurements) {
  ClusterContainer clusters;

  if (!m_cfg.inputClusters.empty()) {
    ACTS_VERBOSE("Fetch clusters for writing: " << m_cfg.inputClusters);
    clusters = m_inputClusters(ctx);
  }

  ACTS_VERBOSE("Writing " << measurements.size()
                          << " measurements in this event.");

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    const auto& from = measurements[hitIdx];
    const Cluster* fromCluster = clusters.empty() ? nullptr : &clusters[hitIdx];

    auto to = m_trackerHitPlaneCollection->create();
    EDM4hepUtil::writeMeasurement(
        from, to, fromCluster, *m_trackerHitRawCollection,
        [](Acts::GeometryIdentifier id) { return id.value(); });
  }

  m_writer.writeEvent();
  m_store.clearCollections();

  return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
