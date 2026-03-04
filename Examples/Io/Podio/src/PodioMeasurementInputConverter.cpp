// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioMeasurementInputConverter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsPodioEdm/MeasurementCollection.h"

#include <podio/Frame.h>

namespace ActsExamples {

PodioMeasurementInputConverter::PodioMeasurementInputConverter(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : PodioInputConverter{"PodioMeasurementInputConverter", config.inputFrame,
                          std::move(logger)},
      m_cfg(config) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument(
        "PodioMeasurementInputConverter: "
        "inputMeasurements is not set");
  }

  m_outputMeasurements.initialize(m_cfg.outputMeasurements);

  m_outputMeasurementParticlesMap.initialize(
      m_cfg.outputMeasurementParticlesMap);
  m_outputParticleMeasurementsMap.initialize(
      m_cfg.outputParticleMeasurementsMap);

  m_outputMeasurementSimHitsMap.initialize(m_cfg.outputMeasurementSimHitsMap);
  m_outputSimHitMeasurementsMap.initialize(m_cfg.outputSimHitMeasurementsMap);

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputSimHitAssociation.initialize(m_cfg.inputSimHitAssociation);
}

ProcessCode PodioMeasurementInputConverter::convert(
    const AlgorithmContext& ctx, const podio::Frame& frame) const {
  const auto& inputMeasurements =
      frame.get<ActsPodioEdm::MeasurementCollection>(m_cfg.inputMeasurements);

  ACTS_DEBUG("Read " << inputMeasurements.size() << " measurements from PODIO");

  MeasurementContainer outputMeasurements;

  const auto& simHits = m_inputSimHits(ctx);
  const auto& simHitAssociations = m_inputSimHitAssociation(ctx);

  IndexMultimap<SimBarcode> measurementToParticles;
  measurementToParticles.reserve(inputMeasurements.size());
  IndexMultimap<Index> measurementToSimHits;
  measurementToParticles.reserve(inputMeasurements.size());

  for (const auto& inputMeas : inputMeasurements) {
    std::size_t index = outputMeasurements.addMeasurement(
        inputMeas.size(), Acts::GeometryIdentifier{inputMeas.getGeometryId()});

    auto meas = outputMeasurements.at(index);

    std::copy(inputMeas.parameterValues_begin(),
              inputMeas.parameterValues_end(), meas.parameters().begin());

    std::span<double> cov{meas.covariance().data(), meas.size() * meas.size()};

    std::copy(inputMeas.covarianceValues_begin(),
              inputMeas.covarianceValues_end(), cov.begin());

    meas.setSubspaceIndices(inputMeas.getIndices());

    // @TODO: Check if this SimHit is naively ordered such that simHit.nth() produces the correct internal sim hit
    // index -> measToSimHit -> simHit -> compare with what we get from edm4hep
    const auto& externalSimHit = inputMeas.getSimHit();

    auto simHitIdx = simHitAssociations.lookup(externalSimHit);
    auto internalSimHitIt = simHits.nth(simHitIdx);
    if (internalSimHitIt == simHits.end()) {
      throw std::invalid_argument("Invalid sim hit index");
    }
    const auto& internalSimHit = *internalSimHitIt;

    ACTS_VERBOSE("Hit lookup ext -> int "
                 << externalSimHit.id() << " -> " << "#" << simHitIdx << " ("
                 << simHitIdx << ") " << internalSimHit.position().transpose());

    measurementToSimHits.emplace_hint(measurementToSimHits.end(), index,
                                      simHitIdx);
    measurementToParticles.emplace_hint(measurementToParticles.end(), index,
                                        internalSimHit.particleId());
  }

  m_outputMeasurements(ctx, std::move(outputMeasurements));

  m_outputParticleMeasurementsMap(ctx,
                                  invertIndexMultimap(measurementToParticles));
  m_outputMeasurementParticlesMap(ctx, std::move(measurementToParticles));

  m_outputSimHitMeasurementsMap(ctx, invertIndexMultimap(measurementToSimHits));
  m_outputMeasurementSimHitsMap(ctx, std::move(measurementToSimHits));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
