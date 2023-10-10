// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <list>
#include <stdexcept>

#include <edm4hep/TrackerHit.h>
#include <edm4hep/TrackerHitCollection.h>
#include <edm4hep/TrackerHitPlane.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepMeasurementReader::EDM4hepMeasurementReader(
    const EDM4hepMeasurementReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepMeasurementReader", level)) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }

  m_reader.openFile(m_cfg.inputPath);

  m_eventsRange = std::make_pair(0, m_reader.getEntries("events"));

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);
  m_outputClusters.maybeInitialize(m_cfg.outputClusters);
}

std::string EDM4hepMeasurementReader::EDM4hepMeasurementReader::name() const {
  return "EDM4hepMeasurementReader";
}

std::pair<size_t, size_t> EDM4hepMeasurementReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode EDM4hepMeasurementReader::read(const AlgorithmContext& ctx) {
  MeasurementContainer measurements;
  ClusterContainer clusters;
  // TODO what about those?
  IndexMultimap<Index> measurementSimHitsMap;
  IndexSourceLinkContainer sourceLinks;

  podio::Frame frame = m_reader.readEntry("events", ctx.eventNumber);

  const auto& trackerHitPlaneCollection =
      frame.get<edm4hep::TrackerHitPlaneCollection>("ActsTrackerHitsPlane");
  const auto& trackerHitRawCollection =
      frame.get<edm4hep::TrackerHitCollection>("ActsTrackerHitsRaw");

  for (const auto& trackerHitPlane : trackerHitPlaneCollection) {
    Cluster cluster;
    auto measurement = EDM4hepUtil::readMeasurement(
        trackerHitPlane, &trackerHitRawCollection, &cluster,
        [](std::uint64_t cellId) { return Acts::GeometryIdentifier(cellId); });

    measurements.push_back(std::move(measurement));
    clusters.push_back(std::move(cluster));
  }

  // Write the data to the EventStore
  m_outputMeasurements(ctx, std::move(measurements));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementSimHitsMap));
  m_outputSourceLinks(ctx, std::move(sourceLinks));
  if (not m_cfg.outputClusters.empty()) {
    m_outputClusters(ctx, std::move(clusters));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
