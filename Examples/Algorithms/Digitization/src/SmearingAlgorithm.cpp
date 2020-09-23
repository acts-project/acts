// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

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
}

ActsExamples::ProcessCode ActsExamples::SmearingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  if (not m_cfg.configured) {
    ACTS_FATAL("Smearing Algorithm is misconfigured. Aborting.");
    return ProcessCode::ABORT;
  }

  // Prepare the input and output collections
  const auto& hits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);

  ActsExamples::GeometryIdMultimap<Acts::FittableMeasurement<DigitizedHit>>
      measurements;

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  for (auto&& [moduleGeoId, moduleHits] : groupByModule(hits)) {
    for (auto ih = moduleHits.begin(); ih != moduleHits.end(); ++ih) {
      const auto& hit = *ih;

      auto surfaceRange = selectModule(m_digitizableSurfaces, moduleGeoId);
      auto smearerRange = selectModule(m_cfg.smearers, moduleGeoId);
      if (not surfaceRange.empty() and not smearerRange.empty()) {
        // First one wins (there shouldn't be more than one smearer per) surface
        auto& surface = surfaceRange.begin()->second;
        auto& smearer = smearerRange.begin()->second;
        ActsFatras::SmearInput sInput(hit, ctx.geoContext, surface);
        // Run the visitor
        std::visit(
            [&](auto&& sm) {
              auto sParSet = sm.first(sInput, rng, sm.second);
              if (sParSet.ok()) {
                auto measurement = createMeasurement(std::move(sParSet.value()),
                                                     surface, {hit});
                measurements.emplace_hint(measurements.end(),
                                          surface->geometryId(),
                                          std::move(measurement));
              }
            },
            smearer);
      }
    }
  }

  // write the clusters to the EventStore
  ctx.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
  return ActsExamples::ProcessCode::SUCCESS;
}