// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsFatras/Digitization/UncorrelatedHitSmearer.hpp"

#include <stdexcept>

ActsExamples::SmearingAlgorithm::SmearingAlgorithm(
    ActsExamples::SmearingAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("SmearingAlgorithm", lvl),
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
  if (not m_cfg.configured) {
    throw std::invalid_argument("Smearing Algorithm is misconfigured");
  }
}

namespace {

/// Create a fittable measurmement from a parameter set
///
/// @param sourceLink The corresponding source link.
/// @param surface The surface of the measurement.
/// @param paramSet The ParameterSet created from the smearer.
/// @return A fittable framework measurement
template <Acts::BoundIndices... kParameters>
ActsExamples::Measurement makeMeasurement(
    ActsExamples::IndexSourceLink sourceLink, const Acts::Surface& surface,
    Acts::ParameterSet<Acts::BoundIndices, kParameters...>&& paramSet) {
  using ConcreteMeasurement =
      Acts::Measurement<ActsExamples::IndexSourceLink, Acts::BoundIndices,
                        kParameters...>;
  return ConcreteMeasurement(surface.getSharedPtr(), sourceLink,
                             std::move(paramSet));
}

}  // namespace

ActsExamples::ProcessCode ActsExamples::SmearingAlgorithm::execute(
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
    const auto& moduleSimHits = simHitsGroup.second;

    auto smearerItr = m_cfg.smearers.find(moduleGeoId);
    if (smearerItr == m_cfg.smearers.end()) {
      ACTS_DEBUG("No smearing function present for module " << moduleGeoId);
      continue;
    }
    const Acts::Surface* surface =
        m_cfg.trackingGeometry->findSurface(moduleGeoId);
    if (not surface) {
      // this is either an invalid geometry id or a misconfigured smearer
      // setup; both cases can not be handled and should be fatal.
      ACTS_ERROR("Could not find surface " << moduleGeoId
                                           << " for configured smearer");
      return ProcessCode::ABORT;
    }

    for (auto ih = moduleSimHits.begin(); ih != moduleSimHits.end(); ++ih) {
      const auto& simHit = *ih;
      const auto simHitIdx = simHits.index_of(ih);

      // run the smearer
      std::visit(
          [&](auto&& smearer) {
            ActsFatras::SmearInput smearInput(simHit, ctx.geoContext, surface);
            auto smearResult = smearer.first(smearInput, rng, smearer.second);
            if (not smearResult.ok()) {
              // do not store un-smearable measurements
              return;
            }

            // the measurement container is unordered and the index under which
            // the measurement will be stored is known before adding it.
            Index hitIdx = measurements.size();
            IndexSourceLink sourceLink(moduleGeoId, hitIdx);
            auto meas = makeMeasurement(sourceLink, *surface,
                                        std::move(smearResult.value()));

            // add to output containers
            // index map and source link container are geometry-ordered.
            // since the input is also geometry-ordered, new items can
            // be added at the end.
            sourceLinks.emplace_hint(sourceLinks.end(), std::move(sourceLink));
            measurements.emplace_back(std::move(meas));
            // this digitization does not do hit merging so there is only one
            // mapping entry for each digitized hit.
            hitParticlesMap.emplace_hint(hitParticlesMap.end(), hitIdx,
                                         simHit.particleId());
            hitSimHitsMap.emplace_hint(hitSimHitsMap.end(), hitIdx, simHitIdx);
          },
          *smearerItr);
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
