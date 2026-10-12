// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/CalorimeterDigitizationAlgorithm.hpp"

#include <stdexcept>
#include <utility>

namespace ActsExamples {

CalorimeterDigitizationAlgorithm::CalorimeterDigitizationAlgorithm(
    const Config& config, Acts::Logging::Level level)
    : IAlgorithm("CalorimeterDigitization",
                 Acts::getDefaultLogger("CalorimeterDigitization", level)),
      m_config(config),
      m_response(config.response) {
  if (config.inputSimHits.empty() || config.outputHits.empty() ||
      config.inputSimHits == config.outputHits) {
    throw std::invalid_argument(
        "Calorimeter input/output names must be distinct "
        "and nonempty");
  }
  m_input.initialize(config.inputSimHits);
  m_output.initialize(config.outputHits);
}

ProcessCode CalorimeterDigitizationAlgorithm::execute(
    const AlgorithmContext& context) const {
  auto hits = m_response(m_input(context));
  ACTS_DEBUG("Produced " << hits.size() << " calibrated calorimeter cells");
  m_output(context, std::move(hits));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
