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

#include "edm4hep/MCParticleCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

namespace ActsExamples {

/// Write particles to EDM4hep
///
/// Inpersistent information:
/// - particle ID
/// - process
class EDM4hepParticleWriter final : public WriterT<SimParticleContainer> {
 public:
  struct Config {
    /// Input particles collection to write.
    std::string inputParticles;
    /// Where to place the output file.
    std::string outputPath;
    /// Name of the particle collection in EDM4hep.
    std::string outputParticles = "MCParticles";
  };

  /// Construct the particle writer.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  EDM4hepParticleWriter(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode endRun() final;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] particles are the particle to be written
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const SimParticleContainer& particles) final;

 private:
  Config m_cfg;

  podio::ROOTWriter m_writer;
  podio::EventStore m_store;

  edm4hep::MCParticleCollection* m_mcParticleCollection;
};

}  // namespace ActsExamples
