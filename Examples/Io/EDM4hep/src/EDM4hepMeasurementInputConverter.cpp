// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementInputConverter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/EDM4hep/TrackerHitCompatibility.hpp"
#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMeasurementInputConverter::EDM4hepMeasurementInputConverter(
    const EDM4hepMeasurementInputConverter::Config& config,
    Acts::Logging::Level level)
    : EDM4hepInputConverter("EDM4hepMeasurementInputConverter", level,
                            config.inputFrame),
      m_cfg(config) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }

  if (m_cfg.inputFrame.empty()) {
    throw std::invalid_argument("Missing input frame");
  }

  m_inputFrame.initialize(m_cfg.inputFrame);

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);
}

ProcessCode EDM4hepMeasurementInputConverter::convert(
    const AlgorithmContext& ctx, const podio::Frame& frame) const {
  MeasurementContainer measurements;
  ClusterContainer clusters;
  // TODO what about those?
  IndexMultimap<Index> measurementSimHitsMap;

  const auto& trackerHitPlaneCollection =
      frame.get<edm4hep::TrackerHitPlaneCollection>("ActsTrackerHitsPlane");
  const auto& trackerHitRawCollection =
      frame.get<edm4hep::TrackerHit3DCollection>("ActsTrackerHitsRaw");

  for (const auto& trackerHitPlane : trackerHitPlaneCollection) {
    Cluster cluster;
    EDM4hepUtil::readMeasurement(
        measurements, trackerHitPlane, &trackerHitRawCollection, &cluster,
        [](std::uint64_t cellId) { return Acts::GeometryIdentifier(cellId); });

    clusters.push_back(std::move(cluster));
  }

  // Write the data to the EventStore
  m_outputMeasurements(ctx, std::move(measurements));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));
  if (!m_cfg.outputClusters.empty()) {
    m_outputClusters(ctx, std::move(clusters));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
