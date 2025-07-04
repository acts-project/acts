// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ParticleSelector.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>

namespace ActsExamples {

bool ParticleSelector::MeasurementCounter::isValidParticle(
    const SimParticle& particle,
    const InverseMultimap<SimBarcode>& particleMeasurementsMap,
    const MeasurementContainer& measurements) const {
  // No hit cuts, accept everything
  if (counters.empty()) {
    return true;
  }

  const auto [measurementsBegin, measurementsEnd] =
      particleMeasurementsMap.equal_range(particle.particleId());

  boost::container::small_vector<unsigned int, 4> counterValues;
  counterValues.resize(counters.size(), 0);

  for (auto measurementIt = measurementsBegin; measurementIt != measurementsEnd;
       ++measurementIt) {
    const auto measurementIndex = measurementIt->second;
    const auto measurement = measurements.at(measurementIndex);

    const auto geoId = measurement.geometryId();

    for (std::size_t i = 0; i < counters.size(); i++) {
      const auto& [counterMap, threshold] = counters[i];
      if (const auto it = counterMap.find(geoId); it != counterMap.end()) {
        counterValues[i]++;
      }
    }
  }

  for (std::size_t i = 0; i < counters.size(); i++) {
    const auto& [counterMap, threshold] = counters[i];
    const unsigned int value = counterValues[i];
    if (value < threshold) {
      return false;
    }
  }

  return true;
}

ParticleSelector::ParticleSelector(const Config& config,
                                   Acts::Logging::Level level)
    : IAlgorithm("ParticleSelector", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.maybeInitialize(
      m_cfg.inputParticleMeasurementsMap);
  m_inputMeasurements.maybeInitialize(m_cfg.inputMeasurements);
  m_outputParticles.initialize(m_cfg.outputParticles);

  if (!m_inputParticleMeasurementsMap.isInitialized() &&
      (m_cfg.measurementsMin > 0 ||
       m_cfg.measurementsMax < std::numeric_limits<std::size_t>::max())) {
    throw std::invalid_argument(
        "Measurement-based cuts require the inputMeasurementParticlesMap");
  }
  if (!m_cfg.measurementCounter.counters.empty() &&
      (!m_inputParticleMeasurementsMap.isInitialized() ||
       !m_inputMeasurements.isInitialized())) {
    throw std::invalid_argument(
        "Measurement count-based cuts require the inputMeasurementParticlesMap "
        "and inputMeasurements");
  }

  ACTS_DEBUG("selection particle rho [" << m_cfg.rhoMin << "," << m_cfg.rhoMax
                                        << ")");
  ACTS_DEBUG("selection particle |z| [" << m_cfg.absZMin << "," << m_cfg.absZMax
                                        << ")");
  ACTS_DEBUG("selection particle time [" << m_cfg.timeMin << ","
                                         << m_cfg.timeMax << ")");
  ACTS_DEBUG("selection particle phi [" << m_cfg.phiMin << "," << m_cfg.phiMax
                                        << ")");
  ACTS_DEBUG("selection particle eta [" << m_cfg.etaMin << "," << m_cfg.etaMax
                                        << ")");
  ACTS_DEBUG("selection particle |eta| [" << m_cfg.absEtaMin << ","
                                          << m_cfg.absEtaMax << ")");
  ACTS_DEBUG("selection particle pt [" << m_cfg.ptMin << "," << m_cfg.ptMax
                                       << ")");
  ACTS_DEBUG("selection particle m [" << m_cfg.mMin << "," << m_cfg.mMax
                                      << ")");
  ACTS_DEBUG("selection particle hits [" << m_cfg.hitsMin << ","
                                         << m_cfg.hitsMax << ")");
  ACTS_DEBUG("selection particle measurements ["
             << m_cfg.measurementsMin << "," << m_cfg.measurementsMax << ")");
  ACTS_DEBUG("remove charged particles " << m_cfg.removeCharged);
  ACTS_DEBUG("remove neutral particles " << m_cfg.removeNeutral);
  ACTS_DEBUG("remove secondary particles " << m_cfg.removeSecondaries);
  ACTS_DEBUG("exclude pdgs: ");
  for (auto pdg : m_cfg.excludeAbsPdgs) {
    ACTS_DEBUG("  " << pdg);
  }
  ACTS_DEBUG("primary vertex ID [" << m_cfg.minPrimaryVertexId << ","
                                   << m_cfg.maxPrimaryVertexId << ")");
}

ProcessCode ParticleSelector::execute(const AlgorithmContext& ctx) const {
  // prepare input/ output types
  const SimParticleContainer& inputParticles = m_inputParticles(ctx);

  const static InverseMultimap<SimBarcode> emptyMeasurementParticlesMap;
  const InverseMultimap<SimBarcode>& inputMeasurementParticlesMap =
      m_inputParticleMeasurementsMap.isInitialized()
          ? m_inputParticleMeasurementsMap(ctx)
          : emptyMeasurementParticlesMap;

  const static MeasurementContainer emptyMeasurements;
  const MeasurementContainer& inputMeasurements =
      m_inputMeasurements.isInitialized() ? m_inputMeasurements(ctx)
                                          : emptyMeasurements;

  std::size_t nInvalidCharge = 0;
  std::size_t nInvalidHitCount = 0;
  std::size_t nInvalidMeasurementCount = 0;
  std::size_t nInvalidMeasurementRegionCount = 0;

  // helper functions to select tracks
  auto within = [](auto x, auto min, auto max) {
    return (min <= x) && (x < max);
  };

  auto isValidParticle = [&](const SimParticle& p) {
    const auto eta = Acts::VectorHelpers::eta(p.direction());
    const auto phi = Acts::VectorHelpers::phi(p.direction());
    const auto rho = Acts::VectorHelpers::perp(p.position());
    // define charge selection
    const bool validNeutral = (p.charge() == 0) && !m_cfg.removeNeutral;
    const bool validCharged = (p.charge() != 0) && !m_cfg.removeCharged;
    const bool validCharge = validNeutral || validCharged;
    const bool validSecondary = !m_cfg.removeSecondaries || !p.isSecondary();
    const bool validPrimaryVertexId =
        within(p.particleId().vertexPrimary(), m_cfg.minPrimaryVertexId,
               m_cfg.maxPrimaryVertexId);

    nInvalidCharge += static_cast<std::size_t>(!validCharge);

    const bool validHitCount =
        within(p.numberOfHits(), m_cfg.hitsMin, m_cfg.hitsMax);
    nInvalidHitCount += static_cast<std::size_t>(!validHitCount);

    const std::size_t measurementCount =
        inputMeasurementParticlesMap.count(p.particleId());
    const bool validMeasurementCount =
        within(measurementCount, m_cfg.measurementsMin, m_cfg.measurementsMax);
    nInvalidMeasurementCount +=
        static_cast<std::size_t>(!validMeasurementCount);

    const bool validMeasurementRegionCount =
        m_cfg.measurementCounter.isValidParticle(
            p, inputMeasurementParticlesMap, inputMeasurements);
    nInvalidMeasurementRegionCount +=
        static_cast<std::size_t>(!validMeasurementRegionCount);

    // Pdg selection
    bool validPdg = true;
    for (auto pdg : m_cfg.excludeAbsPdgs) {
      if (p.absolutePdg() == std::abs(pdg)) {
        validPdg = false;
        break;
      }
    }

    return validPdg && validCharge && validSecondary && validPrimaryVertexId &&
           validHitCount && validMeasurementCount &&
           validMeasurementRegionCount &&
           within(p.transverseMomentum(), m_cfg.ptMin, m_cfg.ptMax) &&
           within(std::abs(eta), m_cfg.absEtaMin, m_cfg.absEtaMax) &&
           within(eta, m_cfg.etaMin, m_cfg.etaMax) &&
           within(phi, m_cfg.phiMin, m_cfg.phiMax) &&
           within(std::abs(p.position()[Acts::ePos2]), m_cfg.absZMin,
                  m_cfg.absZMax) &&
           within(rho, m_cfg.rhoMin, m_cfg.rhoMax) &&
           within(p.time(), m_cfg.timeMin, m_cfg.timeMax) &&
           within(p.mass(), m_cfg.mMin, m_cfg.mMax);
  };

  SimParticleContainer outputParticles;
  outputParticles.reserve(inputParticles.size());

  // copy selected particles
  for (const auto& inputParticle : inputParticles) {
    if (!isValidParticle(inputParticle)) {
      continue;
    }

    outputParticles.insert(outputParticles.end(), inputParticle);
  }
  outputParticles.shrink_to_fit();

  ACTS_DEBUG("event " << ctx.eventNumber << " selected "
                      << outputParticles.size() << " from "
                      << inputParticles.size() << " particles");
  ACTS_DEBUG("filtered out because of charge: " << nInvalidCharge);
  ACTS_DEBUG("filtered out because of hit count: " << nInvalidHitCount);
  ACTS_DEBUG("filtered out because of measurement count: "
             << nInvalidMeasurementCount);

  m_outputParticles(ctx, std::move(outputParticles));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
