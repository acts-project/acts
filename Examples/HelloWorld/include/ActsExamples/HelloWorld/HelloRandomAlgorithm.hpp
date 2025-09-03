// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <array>
#include <memory>
#include <string>

#include "HelloData.hpp"

namespace ActsExamples {

/// An example algorithm that uses the random number generator to generate data.
class HelloRandomAlgorithm : public ActsExamples::IAlgorithm {
 public:
  struct Config {
    std::shared_ptr<ActsExamples::RandomNumbers> randomNumbers = nullptr;
    /// Random distribution parameters.
    std::array<double, 2> gaussParameters = {{0., 1.}};
    std::array<double, 2> uniformParameters = {{0., 1.}};
    std::array<double, 2> gammaParameters = {{0., 1.}};
    int poissonParameter = 40;
    std::size_t drawsPerEvent = 0;
    /// Where to store the generated data in the event store.
    std::string output;
  };

  explicit HelloRandomAlgorithm(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  /// Generate random numbers from various distributions.
  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Access to the config struct
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  WriteDataHandle<HelloDataCollection> m_writeHandle{this, "Output"};
};

}  // namespace ActsExamples
