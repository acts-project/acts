// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementOutputConverter.hpp"

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMeasurementOutputConverter::EDM4hepMeasurementOutputConverter(
    const EDM4hepMeasurementOutputConverter::Config& config,
    Acts::Logging::Level level)
    : PodioOutputConverter("EDM4hepMeasurementOutputConverter", level),
      m_cfg(config) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputClusters.maybeInitialize(m_cfg.inputClusters);
  m_outputTrackerHitsPlane.initialize(m_cfg.outputTrackerHitsPlane);
  m_outputTrackerHitsRaw.initialize(m_cfg.outputTrackerHitsRaw);
}

ProcessCode EDM4hepMeasurementOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  ClusterContainer clusters;

  edm4hep::TrackerHitPlaneCollection hitsPlane;
  edm4hep::TrackerHit3DCollection hits;

  if (!m_cfg.inputClusters.empty()) {
    ACTS_VERBOSE("Fetch clusters for writing: " << m_cfg.inputClusters);
    clusters = m_inputClusters(ctx);
  }

  const auto measurements = m_inputMeasurements(ctx);

  ACTS_VERBOSE("Writing " << measurements.size()
                          << " measurements in this event.");

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    ConstVariableBoundMeasurementProxy from =
        measurements.getMeasurement(hitIdx);
    const Cluster* fromCluster = clusters.empty() ? nullptr : &clusters[hitIdx];

    auto to = hitsPlane.create();
    EDM4hepUtil::writeMeasurement(
        from, to, fromCluster, hits,
        [](Acts::GeometryIdentifier id) { return id.value(); });
  }

  m_outputTrackerHitsPlane(ctx, std::move(hitsPlane));
  m_outputTrackerHitsRaw(ctx, std::move(hits));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepMeasurementOutputConverter::collections()
    const {
  return {m_cfg.outputTrackerHitsPlane, m_cfg.outputTrackerHitsRaw};
}

}  // namespace ActsExamples
