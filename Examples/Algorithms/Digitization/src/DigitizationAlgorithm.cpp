// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <algorithm>
#include <stdexcept>
#include <type_traits>

ActsExamples::DigitizationAlgorithm::DigitizationAlgorithm(
    ActsExamples::DigitizationAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("DigitizationAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (not m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }
}

ActsExamples::ProcessCode ActsExamples::DigitizationAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // retrieve input
  const auto& simHits = ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimHits);

  // prepare output containers
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurements;
  IndexMultimap<ActsFatras::Barcode> hitParticlesMap;
  IndexMultimap<Index> hitSimHitsMap;
  sourceLinks.reserve(simHits.size());
  measurements.reserve(simHits.size());
  hitParticlesMap.reserve(simHits.size());
  hitSimHitsMap.reserve(simHits.size());

  // setup random number generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  for (auto simHitsGroup : groupByModule(simHits)) {
    // manual pair unpacking instead of using
    //   auto [moduleGeoId, moduleSimHits] : ...
    // otherwise clang on macos complains that it is unable to capture the local
    // binding in the lambda used for visiting the smearer below.
    Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
   
    const Acts::Surface* surfacePtr =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
   
    if (not surfacePtr) {
      // this is either an invalid geometry id or a misconfigured smearer
      // setup; both cases can not be handled and should be fatal.
      ACTS_ERROR("Could not find surface " << moduleGeoId
                                           << " for configured smearer");
      return ProcessCode::ABORT;
    }

  }
  
  ctx.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  ctx.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                     std::move(hitParticlesMap));
  ctx.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                     std::move(hitSimHitsMap));
  return ProcessCode::SUCCESS;
}
