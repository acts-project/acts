// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace ActsExamples {

/// Group particles into proto vertices using truth information.
class TruthVertexFinder final : public IAlgorithm {
 public:
  struct Config {
    /// The input truth particles that should be used to create proto vertices.
    std::string inputParticles;
    /// The output proto vertices collection.
    std::string outputProtoVertices;
    /// Exclude secondary particles not originating from the primary vertex.
    bool excludeSecondaries = false;
    /// Build separate proto vertices for the secondary particles.
    bool separateSecondaries = false;
  };

  TruthVertexFinder(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  WriteDataHandle<ProtoVertexContainer> m_outputProtoVertices{
      this, "OutputProtoVertices"};
};

}  // namespace ActsExamples
