// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <cstddef>
#include <limits>
#include <string>

namespace ActsExamples {
struct AlgorithmContext;

class CsvExaTrkXGraphWriter final
    : public WriterT<std::pair<std::vector<int64_t>, std::vector<float>>> {
 public:
  struct Config {
    /// Which simulated (truth) hits collection to use.
    std::string inputGraph;
    /// Where to place output files
    std::string outputDir;
    /// Output filename stem.
    std::string outputStem = "exatrkx-graph";
  };

  /// Construct the cluster writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  CsvExaTrkXGraphWriter(const Config& config, Acts::Logging::Level level);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] simHits are the simhits to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::pair<std::vector<int64_t>, std::vector<float>>&
                         graph) override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
