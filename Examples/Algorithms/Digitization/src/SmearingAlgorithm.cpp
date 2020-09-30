// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <stdexcept>

ActsExamples::SmearingAlgorithm::SmearingAlgorithm(
    ActsExamples::SmearingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SmearingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing input hits collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing output measurement collection");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }

  // fill the digitizables map to allow lookup by geometry id only
  m_cfg.trackingGeometry->visitSurfaces([this](const Acts::Surface* surface) {
    this->m_dSurfaces.insert_or_assign(surface->geometryId(), surface);
  });
}

ActsExamples::ProcessCode ActsExamples::SmearingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  if (not m_cfg.configured) {
    ACTS_FATAL("Smearing Algorithm is misconfigured. Aborting.");
    return ProcessCode::ABORT;
  }

  const auto& hits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  ActsExamples::GeometryIdMultimap<Acts::FittableMeasurement<DigitizedHit>>
      measurements;
  measurements.reserve(hits.size());

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  for (auto&& [moduleGeoId, moduleHits] : groupByModule(hits)) {
    for (auto ih = moduleHits.begin(); ih != moduleHits.end(); ++ih) {
      // Gather the hit index for further storage
      unsigned int hitidx = ih - hits.begin();
      const auto& hit = *ih;

      auto smearItr = m_cfg.smearers.find(moduleGeoId);
      if (smearItr != m_cfg.smearers.end()) {
        auto surfaceItr = m_dSurfaces.find(moduleGeoId);
        if (surfaceItr != m_dSurfaces.end()) {
          // First one wins (there shouldn't be more than one smearer per)
          // surface
          auto& surface = surfaceItr->second;
          auto& smearer = *smearItr;
          ActsFatras::SmearInput sInput(hit, ctx.geoContext, surface);
          // Run the visitor
          std::visit(
              [&](auto&& sm) {
                auto sParSet = sm.first(sInput, rng, sm.second);
                if (sParSet.ok()) {
                  auto measurement = createMeasurement(
                      std::move(sParSet.value()), *surface, {hitidx});
                  measurements.emplace_hint(measurements.end(),
                                            surface->geometryId(),
                                            std::move(measurement));
                }
              },
              smearer);
        } else {
          ACTS_WARNING("Could not find surface with geometry identifier "
                       << moduleGeoId.value());
        }
      } else {
        ACTS_DEBUG("No smearingm function present for this module in volume "
                   << moduleGeoId.volume());
      }
    }
  }

  // write the clusters to the EventStore
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  return ActsExamples::ProcessCode::SUCCESS;
}
