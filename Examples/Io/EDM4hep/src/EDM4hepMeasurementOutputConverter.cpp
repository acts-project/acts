// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepMeasurementOutputConverter.hpp"

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"
#include <ActsPodioEdm/TrackerHitLocalCollection.h>
#include <ActsPodioEdm/TrackerHitLocalSimTrackerHitLinkCollection.h>

#include <stdexcept>

namespace ActsExamples {

EDM4hepMeasurementOutputConverter::EDM4hepMeasurementOutputConverter(
    const EDM4hepMeasurementOutputConverter::Config& config,
    std::unique_ptr<const Acts::Logger> logger)
    : PodioOutputConverter("EDM4hepMeasurementOutputConverter",
                           std::move(logger)),
      m_cfg(config) {
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::runtime_error(
        "EDM4hepMeasurementOutputConverter: trackingGeometry is null");
  }

  const bool hasAssoc = m_cfg.inputSimHitAssociation.has_value();
  const bool hasMap = m_cfg.inputMeasurementSimHitsMap.has_value();
  const bool hasLinks = m_cfg.outputSimHitLinks.has_value();
  if (hasAssoc != hasMap || hasMap != hasLinks) {
    throw std::invalid_argument(
        "EDM4hepMeasurementOutputConverter: inputSimHitAssociation, "
        "inputMeasurementSimHitsMap, and outputSimHitLinks must all be set or "
        "all be unset");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputTrackerHitsLocal.initialize(m_cfg.outputTrackerHitsLocal);
  m_inputSimHitAssociation.maybeInitialize(m_cfg.inputSimHitAssociation);
  m_inputMeasurementSimHitsMap.maybeInitialize(
      m_cfg.inputMeasurementSimHitsMap);
  m_outputSimHitLinks.maybeInitialize(m_cfg.outputSimHitLinks);
}

ProcessCode EDM4hepMeasurementOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const auto& measurements = m_inputMeasurements(ctx);

  ACTS_VERBOSE("Writing " << measurements.size()
                          << " measurements in this event.");

  auto hits = std::make_unique<ActsPodioEdm::TrackerHitLocalCollection>();

  for (Index hitIdx = 0u; hitIdx < measurements.size(); ++hitIdx) {
    ConstVariableBoundMeasurementProxy from =
        measurements.getMeasurement(hitIdx);

    const Acts::Surface* surface =
        m_cfg.trackingGeometry->findSurface(from.geometryId());
    if (surface == nullptr) {
      throw std::runtime_error(
          "EDM4hepMeasurementOutputConverter: surface not found for geometry "
          "id " +
          std::to_string(from.geometryId().value()));
    }

    auto to = hits->create();
    EDM4hepUtil::writeMeasurement(ctx.geoContext, from, to, *surface);
  }

  if (m_outputSimHitLinks.isInitialized()) {
    const auto& simHitAssociation = m_inputSimHitAssociation(ctx);
    const auto& measToSimHits = m_inputMeasurementSimHitsMap(ctx);

    ActsPlugins::EDM4hepUtil::SimHitForHitIndex lookup =
        [&measToSimHits,
         &simHitAssociation](std::size_t i) -> std::optional<edm4hep::SimTrackerHit> {
      auto it = measToSimHits.find(static_cast<Index>(i));
      if (it == measToSimHits.end()) {
        return std::nullopt;
      }
      return simHitAssociation.lookup(it->second);
    };

    auto links = std::make_unique<
        ActsPodioEdm::TrackerHitLocalSimTrackerHitLinkCollection>();
    ActsPlugins::EDM4hepUtil::writeTrackerHitSimHitLinks(*hits, *links, lookup);
    m_outputSimHitLinks(ctx, std::move(links));
  }

  m_outputTrackerHitsLocal(ctx, std::move(hits));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepMeasurementOutputConverter::collections()
    const {
  std::vector<std::string> result = {m_cfg.outputTrackerHitsLocal};
  if (m_cfg.outputSimHitLinks.has_value()) {
    result.push_back(m_cfg.outputSimHitLinks.value());
  }
  return result;
}

}  // namespace ActsExamples
