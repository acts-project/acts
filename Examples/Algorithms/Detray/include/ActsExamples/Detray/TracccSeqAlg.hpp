// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>

// Forward-declare to keep CUDA/traccc headers out of this header.
namespace ActsExamples {
struct TracccChain;
}

namespace ActsExamples {

class TracccSeqAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Path to the detector geometry file
    std::string detectorFile;
    /// Path to the digitization config file
    std::string digitizationFile;
    /// Path to the detector conditions file (optional)
    std::string conditionsFile;
    /// Path to the material file (optional)
    std::string materialFile;
    /// Path to the grid file (optional)
    std::string gridFile;
    /// Path to the magnetic field file
    std::string bfieldFile;
    /// Directory containing per-event cell CSV files
    std::string dataDirectory;
  };

  explicit TracccSeqAlgorithm(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const override;
  ProcessCode finalize() override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::shared_ptr<TracccChain> m_chain;
};

}  // namespace ActsExamples
