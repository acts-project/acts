// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioMeasurementOutputConverter.hpp"

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"
#include "ActsPodioEdm/MeasurementCollection.h"

#include <cstdint>

#include <podio/CollectionBase.h>

namespace ActsExamples {

PodioMeasurementOutputConverter::PodioMeasurementOutputConverter(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : PodioOutputConverter{"PodioMeasurementOutputConverter",
                           std::move(logger)},
      m_cfg{config} {
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);

  if (m_cfg.inputSimHitAssociation.has_value() !=
      m_cfg.inputMeasurementSimHitsMap.has_value()) {
    throw std::invalid_argument(
        "PodioMeasurementOutputConverter: inputSimHitAssociation and "
        "inputMeasurementSimHitsMap must be both set or both unset");
  }

  m_inputSimHitAssociation.maybeInitialize(m_cfg.inputSimHitAssociation);
  m_inputMeasurementSimHitsMap.maybeInitialize(
      m_cfg.inputMeasurementSimHitsMap);
}

ProcessCode PodioMeasurementOutputConverter::execute(
    const AlgorithmContext& context) const {
  const auto& measurements = m_inputMeasurements(context);

  const auto& measToSimHits = m_inputMeasurementSimHitsMap(context);
  // Map internal sim hits to edm4hep input sim hits
  const auto& simHitAssociation = m_inputSimHitAssociation(context);

  auto outputMeasurements =
      std::make_unique<ActsPodioEdm::MeasurementCollection>();

  for (const MeasurementContainer::ConstVariableProxy& meas : measurements) {
    auto to = outputMeasurements->create();
    for (const auto val : meas.subspaceIndexVector()) {
      to.addToIndices(static_cast<std::uint16_t>(val));
    }

    for (const double val : meas.parameters()) {
      to.addToParameterValues(val);
    }

    std::span<const double> cov{meas.covariance().data(),
                                meas.size() * meas.size()};

    for (const double val : cov) {
      to.addToCovarianceValues(val);
    }

    to.setGeometryId(meas.geometryId().value());

    throw_assert(to.size() == meas.size(), "Invalid sizes after filling");

    // Get the sim hit associated with this measurement
    const auto simHitIt = measToSimHits.find(meas.index());
    if (simHitIt == measToSimHits.end()) {
      ACTS_ERROR("No sim hit found for measurement index " << meas.index());
      return ProcessCode::ABORT;
    }

    // Get the edm4hep input sim hit from the internal index
    const auto& sourceHit = simHitAssociation.lookup(simHitIt->second);
    to.setSimHit(sourceHit);
  }

  m_outputMeasurements(context, std::move(outputMeasurements));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> PodioMeasurementOutputConverter::collections() const {
  return {m_cfg.outputMeasurements};
}

}  // namespace ActsExamples
