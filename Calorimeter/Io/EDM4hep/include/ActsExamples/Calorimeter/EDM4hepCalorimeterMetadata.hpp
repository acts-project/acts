// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Calorimeter/CalorimeterClusteringAlgorithm.hpp"
#include "ActsExamples/Calorimeter/CalorimeterDigitizationAlgorithm.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterClusterOutputConverter.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterInputConverter.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterOutputConverter.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <map>
#include <string>
#include <vector>

#include <podio/Frame.h>

namespace ActsExamples {
/// Attach immutable reconstruction configuration to each event's original
/// frame. Runs after converters; transfers frame ownership to a distinct
/// whiteboard key.
class EDM4hepCalorimeterMetadata final : public IAlgorithm {
 public:
  struct Parameters {
    std::map<std::string, std::vector<std::string>> strings;
    std::map<std::string, std::vector<double>> numbers;
  };
  struct Config {
    std::string inputFrame = "events";
    std::string outputFrame;
    /// Versioned geometry source, including the cell-centre lookup convention.
    std::string geometryIdentifier;
    EDM4hepCalorimeterInputConverter::Config input;
    CalorimeterDigitizationAlgorithm::Config response;
    CalorimeterClusteringAlgorithm::Config clustering;
    EDM4hepCalorimeterOutputConverter::Config hitOutput;
    EDM4hepCalorimeterClusterOutputConverter::Config clusterOutput;
    /// Jet converter metadata, when jets are enabled. Keys must be unique and
    /// start with acts.calo.; existing event keys under that prefix are
    /// rejected.
    Parameters additional;
  };
  explicit EDM4hepCalorimeterMetadata(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);
  ProcessCode execute(const AlgorithmContext& context) const override;

 private:
  Parameters m_parameters;
  ConsumeDataHandle<podio::Frame> m_input{this, "InputFrame"};
  WriteDataHandle<podio::Frame> m_output{this, "OutputFrame"};
};
}  // namespace ActsExamples
