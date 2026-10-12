// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/Digitization/CalorimeterResponse.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

/// ActsExamples adapter for the deterministic prototype calorimeter response.
class CalorimeterDigitizationAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputSimHits;
    std::string outputHits;
    ActsCalorimeter::CalorimeterResponse::Config response;
  };

  explicit CalorimeterDigitizationAlgorithm(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  ProcessCode execute(const AlgorithmContext& context) const override;
  const Config& config() const { return m_config; }

 private:
  Config m_config;
  ActsCalorimeter::CalorimeterResponse m_response;
  ReadDataHandle<ActsCalorimeter::SimCalorimeterHitContainer> m_input{
      this, "InputSimHits"};
  WriteDataHandle<ActsCalorimeter::CalorimeterHitContainer> m_output{
      this, "OutputHits"};
};

}  // namespace ActsExamples
