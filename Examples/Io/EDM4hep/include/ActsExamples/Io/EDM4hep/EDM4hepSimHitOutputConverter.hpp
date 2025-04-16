// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepOutputConverter.hpp"
#include "ActsExamples/Io/Podio/CollectionBaseWriteHandle.hpp"

#include <string>

namespace ActsExamples {

/// Write out a simhit collection to EDM4hep objects
///
/// Inpersistent information:
/// - particle ID
/// - after4 momentum
/// - hit index
/// - digitization channel
class EDM4hepSimHitOutputConverter final : public EDM4hepOutputConverter {
 public:
  struct Config {
    /// Which simulated (truth) hits collection to use.
    std::string inputSimHits;
    /// Which simulated (truth) particle collection to use.
    std::string inputParticles;
    /// WWhere to write the output file to.
    /// Name of the particle collection in EDM4hep.
    std::string outputParticles;
    /// Name of the sim tracker hit collection in EDM4hep
    std::string outputSimTrackerHits = "ActsSimTrackerHits";
  };

  /// Construct the cluster writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepSimHitOutputConverter(const Config& config,
                               Acts::Logging::Level level);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Access to the collections
  std::vector<std::string> collections() const override;

  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] simHits are the simhits to be written
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  Config m_cfg;

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  CollectionBaseWriteHandle m_outputParticles{this, "OutputParticles"};
  CollectionBaseWriteHandle m_outputSimTrackerHits{this,
                                                   "OutputSimTrackerHits"};
};

}  // namespace ActsExamples
