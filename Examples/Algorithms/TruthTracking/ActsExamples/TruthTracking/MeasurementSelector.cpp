// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/MeasurementSelector.hpp"

#include "ActsExamples/EventData/Measurement.hpp"

ActsExamples::MeasurementSelector::MeasurementSelector(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("MeasurementSelector", level), m_cfg(config) {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputParticlesSelected.maybeInitialize(m_cfg.inputParticlesSelected);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);

  ACTS_DEBUG("selection particles " << m_cfg.inputParticlesSelected);
  ACTS_DEBUG("selection primary vertex ID [" << m_cfg.minPrimaryVertexId << ","
                                             << m_cfg.maxPrimaryVertexId
                                             << ")");
}

ActsExamples::ProcessCode ActsExamples::MeasurementSelector::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  const MeasurementContainer& measurements = m_inputMeasurements(ctx);
  const auto& measurementParticlesMap = m_inputMeasurementParticlesMap(ctx);
  const SimParticleContainer* particlesSelected =
      m_inputParticlesSelected.isInitialized() ? &m_inputParticlesSelected(ctx)
                                               : nullptr;

  MeasurementContainer selectedMeasurements;
  IndexMultimap<SimBarcode> selectedMeasurementParticlesMap;

  for (const auto& measurement : measurements) {
    const auto& particleIdRange =
        measurementParticlesMap.equal_range(measurement.index());

    bool particleSelected = false;
    bool primaryVertexSelected = false;

    for (auto particleIdIt = particleIdRange.first;
         particleIdIt != particleIdRange.second; ++particleIdIt) {
      const auto& particleId = particleIdIt->second;

      particleSelected |= (particlesSelected == nullptr) ||
                          particlesSelected->contains(particleId);
      primaryVertexSelected |=
          (particleId.vertexPrimary() >= m_cfg.minPrimaryVertexId) &&
          (particleId.vertexPrimary() < m_cfg.maxPrimaryVertexId);
    }

    const bool selected = particleSelected && primaryVertexSelected;

    if (selected) {
      const auto measurementIndex = selectedMeasurements.size();
      selectedMeasurements.copyMeasurement(measurement);

      for (auto particleIdIt = particleIdRange.first;
           particleIdIt != particleIdRange.second; ++particleIdIt) {
        const auto& particleId = particleIdIt->second;

        selectedMeasurementParticlesMap.emplace_hint(
            measurementParticlesMap.end(), measurementIndex, particleId);
      }
    }
  }

  ACTS_DEBUG("selected " << selectedMeasurements.size() << " from "
                         << measurements.size() << " hits");

  m_outputMeasurements(ctx, std::move(selectedMeasurements));
  m_outputMeasurementParticlesMap(ctx,
                                  std::move(selectedMeasurementParticlesMap));

  return ProcessCode::SUCCESS;
}
