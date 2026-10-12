// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/CalorimeterJetAlgorithm.hpp"

#include <stdexcept>
#include <utility>

namespace ActsExamples {

CalorimeterJetAlgorithm::CalorimeterJetAlgorithm(const Config& config,
                                                 Acts::Logging::Level level)
    : IAlgorithm("CalorimeterJets",
                 Acts::getDefaultLogger("CalorimeterJets", level)),
      m_config(config),
      m_reconstruction(config.jets) {
  if (config.inputClusters.empty() || config.outputJets.empty() ||
      config.inputClusters == config.outputJets) {
    throw std::invalid_argument(
        "Calorimeter input/output names must be distinct and nonempty");
  }
  m_input.initialize(config.inputClusters);
  m_output.initialize(config.outputJets);
}

ProcessCode CalorimeterJetAlgorithm::execute(
    const AlgorithmContext& context) const {
  auto jets = m_reconstruction(m_input(context));
  ACTS_DEBUG("Produced " << jets.size() << " calorimeter jets");
  m_output(context, std::move(jets));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
