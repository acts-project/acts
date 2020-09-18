// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace ActsExamples {

/// @brief Digitization Algorithm that turns simulated
/// hits into measuremetns for Fitting
class DigitizationAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits
    std::string inputSimulatedHits;
    /// Output collection of measuremetns
    std::string outputMeasurements;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers;
  };

  /// Construct the digitization algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  DigitizationAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build measurement from simulation hits at input
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  struct Digitizable {};

  Config m_cfg;
  /// Lookup container for all digitizable surfaces
  std::unordered_map<Acts::GeometryIdentifier, Digitizable> m_digitizables;
};

}  // namespace ActsExamples