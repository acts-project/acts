// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/HitSmearing.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

ActsExamples::HitSmearing::HitSmearing(const Config& cfg,
                                       Acts::Logging::Level lvl)
    : BareAlgorithm("HitSmearing", lvl), m_cfg(cfg) {
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits input collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if ((m_cfg.sigmaLoc0 < 0) or (m_cfg.sigmaLoc1 < 0)) {
    throw std::invalid_argument("Invalid resolution setting");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (not m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }
}

ActsExamples::ProcessCode ActsExamples::HitSmearing::execute(
    const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

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
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  // setup local covariance
  Acts::SymMatrix2 cov = Acts::SymMatrix2::Zero();
  cov(0, 0) = m_cfg.sigmaLoc0 * m_cfg.sigmaLoc0;
  cov(1, 1) = m_cfg.sigmaLoc1 * m_cfg.sigmaLoc1;

  for (auto&& [moduleGeoId, moduleSimHits] : groupByModule(simHits)) {
    // check if we should create hits for this surface
    const Acts::Surface* surface =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
    if (not surface) {
      continue;
    }

    // use iterators manually so we can retrieve the hit index in the container
    for (auto ih = moduleSimHits.begin(); ih != moduleSimHits.end(); ++ih) {
      const auto& simHit = *ih;
      const auto simHitIdx = simHits.index_of(ih);

      // transform global position into local coordinates
      auto lpResult = surface->globalToLocal(ctx.geoContext, simHit.position(),
                                             simHit.unitDirection(), 0.5_um);
      if (not lpResult.ok()) {
        ACTS_ERROR("Global to local transformation did not succeed.");
        return ProcessCode::ABORT;
      }

      // create smeared local measurement
      Acts::Vector2 loc = lpResult.value();
      loc[0] += m_cfg.sigmaLoc0 * stdNormal(rng);
      loc[1] += m_cfg.sigmaLoc1 * stdNormal(rng);

      // the measurement container is unordered and the index under which the
      // measurement will be stored is known before adding it.
      Index hitIdx = measurements.size();
      IndexSourceLink sourceLink(moduleGeoId, hitIdx);
      auto meas = Acts::makeMeasurement(sourceLink, loc, cov, Acts::eBoundLoc0,
                                        Acts::eBoundLoc1);

      // add to output containers. since the input is already geometry-order,
      // new elements in geometry containers can just be appended at the end.
      sourceLinks.emplace_hint(sourceLinks.end(), std::move(sourceLink));
      measurements.emplace_back(std::move(meas));
      // no hit merging -> only one mapping per digitized hit.
      hitParticlesMap.emplace_hint(hitParticlesMap.end(), hitIdx,
                                   simHit.particleId());
      hitSimHitsMap.emplace_hint(hitSimHitsMap.end(), hitIdx, simHitIdx);
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
