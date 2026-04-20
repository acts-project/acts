// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/MeasurementFilterAlgorithm.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <stdexcept>

namespace ActsExamples {

MeasurementFilterAlgorithm::MeasurementFilterAlgorithm(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("MeasurementFilterAlgorithm", std::move(logger)),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurements");
  }
  if (m_cfg.inputMeasurementMap.empty()) {
    throw std::invalid_argument("Missing input measurement map");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing output measurements");
  }
  if (m_cfg.outputIndexRemapping.empty()) {
    throw std::invalid_argument("Missing output index remapping");
  }

  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputMeasurementMap.initialize(m_cfg.inputMeasurementMap);
  m_outputMeasurements.initialize(m_cfg.outputMeasurements);
  m_outputIndexRemapping.initialize(m_cfg.outputIndexRemapping);
}

ProcessCode MeasurementFilterAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& inputMeasurements = m_inputMeasurements(ctx);
  const auto& measurementMap = m_inputMeasurementMap(ctx);

  MeasurementContainer filtered;
  filtered.reserve(inputMeasurements.size());

  MeasurementIndexRemapping remapping;
  remapping.reserve(inputMeasurements.size());

  std::size_t nFiltered = 0;

  for (Index originalIdx = 0;
       originalIdx < static_cast<Index>(inputMeasurements.size());
       ++originalIdx) {
    if (measurementMap.count(originalIdx) != 0u) {
      ++nFiltered;
      continue;
    }
    const auto src = inputMeasurements.getMeasurement(originalIdx);
    filtered.copyMeasurement(src);
    remapping.push_back(originalIdx);
  }

  ACTS_DEBUG("Filtered " << nFiltered << " used measurement(s) out of "
                         << inputMeasurements.size() << "; "
                         << filtered.size() << " remaining");

  m_outputMeasurements(ctx, std::move(filtered));
  m_outputIndexRemapping(ctx, std::move(remapping));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
