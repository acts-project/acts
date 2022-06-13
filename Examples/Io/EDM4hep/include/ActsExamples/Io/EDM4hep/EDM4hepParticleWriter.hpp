// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <string>

namespace ActsExamples {

class EDM4hepParticleWriter final : public WriterT<SimParticleContainer> {
 public:
  struct Config {
    /// Input particles collection to write.
    std::string inputParticles;
    /// Where to place the output file.
    std::string outputPath;
  };

  /// Construct the particle writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  EDM4hepParticleWriter(const Config& cfg, Acts::Logging::Level lvl);

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const SimParticleContainer& particles) final override;

 private:
  Config m_cfg;  //!< Nested configuration struct
};

}  // namespace ActsExamples
