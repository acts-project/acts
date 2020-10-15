// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/HitSmearing.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSourceLink.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

using namespace Acts::UnitLiterals;
ActsExamples::HitSmearing::HitSmearing(const Config& cfg,
                                       Acts::Logging::Level lvl)
    : BareAlgorithm("HitSmearing", lvl), m_cfg(cfg) {
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing output source links collection");
  }
  if ((m_cfg.sigmaLoc0 < 0) or (m_cfg.sigmaLoc1 < 0)) {
    throw std::invalid_argument("Invalid resolution setting");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }
}

ActsExamples::ProcessCode ActsExamples::HitSmearing::execute(
    const AlgorithmContext& ctx) const {
  // setup input and output containers
  const auto& hits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);
  SimSourceLinkContainer sourceLinks;
  sourceLinks.reserve(hits.size());

  // setup random number generator
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  // setup local covariance
  // TODO add support for per volume/layer/module settings
  Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
  cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = m_cfg.sigmaLoc0 * m_cfg.sigmaLoc0;
  cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = m_cfg.sigmaLoc1 * m_cfg.sigmaLoc1;

  for (auto&& [moduleGeoId, moduleHits] : groupByModule(hits)) {
    // check if we should create hits for this surface
    const Acts::Surface* surface =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
    if (not surface) {
      continue;
    }

    // smear all truth hits for this module
    for (const auto& hit : moduleHits) {
      // transform global position into local coordinates
      auto lpResult = surface->globalToLocal(ctx.geoContext, hit.position(),
                                             hit.unitDirection(), 0.5_um);
      Acts::Vector2D lp{0., 0.};
      if (not lpResult.ok()) {
        ACTS_ERROR("Global to local transformation did not succeed.");
        return ProcessCode::ABORT;
      } else {
        lp = lpResult.value();
      }

      // smear truth to create local measurement
      Acts::BoundVector loc = Acts::BoundVector::Zero();
      loc[Acts::eBoundLoc0] = lp[0] + m_cfg.sigmaLoc0 * stdNormal(rng);
      loc[Acts::eBoundLoc1] = lp[1] + m_cfg.sigmaLoc1 * stdNormal(rng);

      // create source link at the end of the container
      auto it = sourceLinks.emplace_hint(sourceLinks.end(), *surface, hit, 2,
                                         loc, cov);
      // ensure hits and links share the same order to prevent ugly surprises
      if (std::next(it) != sourceLinks.end()) {
        ACTS_FATAL("The hit ordering broke. Run for your life.");
        return ProcessCode::ABORT;
      }
    }
  }

  ctx.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
  return ProcessCode::SUCCESS;
}
