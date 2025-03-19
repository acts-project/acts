// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"

#include <mutex>
#include <string>

namespace ActsExamples {

/// Write particles to EDM4hep
///
/// Inpersistent information:
/// - particle ID
/// - process
class EDM4hepParticleOutputConverter final : public IAlgorithm {
 public:
  struct Config {
    /// Input particles collection to write.
    std::string inputParticles;
    /// Name of the particle collection in EDM4hep.
    std::string outputParticles = "MCParticles";
  };

  /// Construct the particle writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  EDM4hepParticleOutputConverter(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const ActsExamples::AlgorithmContext& ctx) const final;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const;

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  CollectionBaseWriteHandle m_outputParticles{this, "OutputParticles"};
};

}  // namespace ActsExamples
