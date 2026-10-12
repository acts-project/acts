// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/Jets/CalorimeterJetReconstruction.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace ActsExamples {

/// ActsExamples adapter producing event-owned jets from calibrated clusters.
class CalorimeterJetAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputClusters;
    std::string outputJets;
    ActsCalorimeter::CalorimeterJetReconstruction::Config jets;
  };

  explicit CalorimeterJetAlgorithm(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  ProcessCode execute(const AlgorithmContext& context) const override;
  const Config& config() const { return m_config; }

 private:
  Config m_config;
  ActsCalorimeter::CalorimeterJetReconstruction m_reconstruction;
  ReadDataHandle<ActsCalorimeter::CalorimeterClusterContainer> m_input{
      this, "InputClusters"};
  WriteDataHandle<ActsCalorimeter::CalorimeterJetContainer> m_output{
      this, "OutputJets"};
};

}  // namespace ActsExamples
