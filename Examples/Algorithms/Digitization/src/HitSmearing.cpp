// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Digitization/HitSmearing.hpp"

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"

FW::HitSmearing::HitSmearing(const Config& cfg, Acts::Logging::Level lvl)
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
  // fill the surface map to allow lookup by geometry id only
  m_cfg.trackingGeometry->visitSurfaces([this](const Acts::Surface* surface) {
    // for now we just require a valid surface
    if (not surface) {
      return;
    }
    this->m_surfaces.insert_or_assign(surface->geoID(), surface);
  });
}

FW::ProcessCode FW::HitSmearing::execute(const AlgorithmContext& ctx) const {
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
  cov(Acts::eLOC_0, Acts::eLOC_0) = m_cfg.sigmaLoc0 * m_cfg.sigmaLoc0;
  cov(Acts::eLOC_1, Acts::eLOC_1) = m_cfg.sigmaLoc1 * m_cfg.sigmaLoc1;

  for (auto&& [moduleGeoId, moduleHits] : groupByModule(hits)) {
    // check if we should create hits for this surface
    const auto is = m_surfaces.find(moduleGeoId);
    if (is == m_surfaces.end()) {
      continue;
    }

    // smear all truth hits for this module
    const Acts::Surface* surface = is->second;
    for (const auto& hit : moduleHits) {
      // transform global position into local coordinates
      Acts::Vector2D pos(0, 0);
      surface->globalToLocal(ctx.geoContext, hit.position(),
                             hit.unitDirection(), pos);

      // smear truth to create local measurement
      Acts::BoundVector loc = Acts::BoundVector::Zero();
      loc[Acts::eLOC_0] = pos[0] + m_cfg.sigmaLoc0 * stdNormal(rng);
      loc[Acts::eLOC_1] = pos[1] + m_cfg.sigmaLoc1 * stdNormal(rng);

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
