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

#include <list>

#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

namespace ActsExamples {

EDM4hepMeasurementReader::EDM4hepMeasurementReader(
    const EDM4hepMeasurementReader::Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepMeasurementReader", level)) {
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurement output collection");
  }

  m_reader.openFile(m_cfg.inputPath);
  m_store.setReader(&m_reader);

  m_eventsRange = std::make_pair(0, m_reader.getEntries());

  m_trackerHitPlaneCollection =
      &m_store.get<edm4hep::TrackerHitPlaneCollection>("ActsTrackerHitsPlane");
  m_trackerHitRawCollection =
      &m_store.create<edm4hep::TrackerHitCollection>("ActsTrackerHitsRaw");
}

std::string EDM4hepMeasurementReader::EDM4hepMeasurementReader::name() const {
  return "EDM4hepMeasurementReader";
}

std::pair<size_t, size_t> EDM4hepMeasurementReader::availableEvents() const {
  return m_eventsRange;
}

ProcessCode EDM4hepMeasurementReader::read(const AlgorithmContext& ctx) {
  // Prepare containers for the hit data using the framework event data types
  GeometryIdMultimap<Measurement> orderedMeasurements;
  ClusterContainer clusters;
  IndexMultimap<Index> measurementSimHitsMap;
  IndexSourceLinkContainer sourceLinks;
  // need list here for stable addresses
  std::list<IndexSourceLink> sourceLinkStorage;

  // TODO reserve space for the containers
  // orderedMeasurements.reserve(measurementData.size());
  // Safe long as we have single particle to sim hit association
  // measurementSimHitsMap.reserve(measurementData.size());
  // sourceLinks.reserve(measurementData.size());

  m_store.clear();
  m_reader.goToEvent(ctx.eventNumber);

  // TODO
  // for (const auto& trackerHitPlane : *m_trackerHitPlaneCollection) { }

  MeasurementContainer measurements;
  for (auto& [_, meas] : orderedMeasurements) {
    measurements.emplace_back(std::move(meas));
  }

  // Write the data to the EventStore
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                     std::move(measurementSimHitsMap));
  ctx.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
  ctx.eventStore.add(m_cfg.outputSourceLinks + "__storage",
                     std::move(sourceLinkStorage));
  if (not m_cfg.outputClusters.empty()) {
    ctx.eventStore.add(m_cfg.outputClusters, std::move(clusters));
  }

  m_reader.endOfEvent();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
