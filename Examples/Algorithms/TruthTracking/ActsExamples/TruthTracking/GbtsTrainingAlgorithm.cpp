// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/GbtsTrainingAlgorithm.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

GbtsTrainingAlgorithm::GbtsTrainingAlgorithm(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("GbtsTrainingAlgorithm", std::move(logger)), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }

  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurements collection");
  }
  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing input simulated hits measurements map");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);

  m_gbtsTrainingTool.emplace(
      m_cfg.gbtsTrainingConfig, config.geometryFileDir,
      this->logger().cloneWithSuffix("GbtsLayerConnectionTool"));
}

GbtsTrainingAlgorithm::~GbtsTrainingAlgorithm() {
  m_gbtsTrainingTool->createConnectionTable(m_cfg.outputFileDir);
}
ProcessCode GbtsTrainingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
  const auto& measurementsIn = m_inputMeasurements(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  ACTS_VERBOSE("analysis hit information for " << particles.size()
                                               << " particles");

  for (const auto& [i, particle] : Acts::enumerate(particles)) {
    // find the corresponding hits for this particle
    const auto& measurements =
        makeRange(particleMeasurementsMap.equal_range(particle.particleId()));
    ACTS_VERBOSE(measurements.size()
                 << " measurements for particle " << particle);

    std::vector<double> hitTimes;
    std::vector<Acts::Experimental::GbtsLayerConnectionTool::TrackCoordinates>
        hitCoords;

    hitTimes.reserve(measurements.size());
    hitCoords.reserve(measurements.size());

    for (const auto& [barcode, index] : measurements) {
      ConstVariableBoundMeasurementProxy measurement =
          measurementsIn.getMeasurement(index);

      ACTS_VERBOSE("   - Measurement " << index << " with barcode " << barcode
                                       << " at " << measurement.geometryId());

      const auto simHitMapIt = measurementSimHitsMap.find(index);
      if (simHitMapIt == measurementSimHitsMap.end()) {
        ACTS_WARNING("No sim hit found for measurement index " << index);
        continue;
      }

      const auto simHitIndex = simHitMapIt->second;

      const auto simHitIt = simHits.nth(simHitIndex);
      if (simHitIt == simHits.end()) {
        ACTS_WARNING("No sim hit found for sim hit index "
                     << simHitIndex << " from measurement " << index);
        continue;
      }

      const auto& simHit = *simHitIt;
      const auto pos = simHit.position();

      const float r = static_cast<float>(std::hypot(pos.x(), pos.y()));
      const float z = static_cast<float>(pos.z());

      hitTimes.emplace_back(simHit.time());
      hitCoords.push_back({r, z});
    }

    std::vector<std::size_t> indices;
    indices.resize(hitCoords.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::ranges::sort(indices, [&hitTimes](std::size_t a, std::size_t b) {
      return hitTimes[a] < hitTimes[b];
    });

    std::vector<Acts::Experimental::GbtsLayerConnectionTool::TrackCoordinates>
        coords;

    coords.reserve(hitCoords.size());

    for (const auto& idx : indices) {
      coords.push_back(hitCoords[idx]);
    }

    {
      std::lock_guard<std::mutex> lock(m_gbtsTrainingToolMutex);
      m_gbtsTrainingTool->addTrack(coords);
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
