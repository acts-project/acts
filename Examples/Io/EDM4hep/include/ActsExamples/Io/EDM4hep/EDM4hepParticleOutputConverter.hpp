// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepOutputConverter.hpp"
#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"

#include <string>

namespace ActsExamples {

/// Write particles to EDM4hep objects
///
/// Inpersistent information:
/// - particle ID
/// - process
class EDM4hepParticleOutputConverter final : public EDM4hepOutputConverter {
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

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  std::vector<std::string> collections() const override;

 private:
  ProcessCode execute(const ActsExamples::AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  CollectionBaseWriteHandle m_outputParticles{this, "OutputParticles"};
};

}  // namespace ActsExamples
