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
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Traccc/TracccChain.hpp"

#include <string>

#include "traccc/edm/silicon_cell_collection.hpp"

namespace ActsExamples::Traccc {

class TracccSeqAlg final : public IAlgorithm {
 public:
  struct Config {
    /// Full traccc chain configuration (geometry, seeding params, etc.)
    TracccChainConfig chain;

    /// WhiteBoard key for the input collections
    std::string inputMeasurements = "";
    std::string inputSpacepoints = "";
  };

  TracccSeqAlg(Config cfg, Acts::Logging::Level logLevel);
  ~TracccSeqAlg() override = default;

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  TracccChain m_chain;

  ReadDataHandle<traccc::edm::measurement_collection::host> m_inputMeasurements{
      this, "inputMeasurements"};
  ReadDataHandle<traccc::edm::spacepoint_collection::host> m_inputSpacepoints{
      this, "inputSpacepoints"};
};

}  // namespace ActsExamples::Traccc
