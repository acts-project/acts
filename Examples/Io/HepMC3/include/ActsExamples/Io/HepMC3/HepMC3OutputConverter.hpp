// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace HepMC3 {
class GenEvent;
}

namespace ActsExamples {

class HepMC3OutputConverter : public IAlgorithm {
 public:
  struct Config {
    std::string inputParticles;
    std::string inputVertices;
    std::string outputEvent;
  };
  HepMC3OutputConverter(const Config& config, Acts::Logging::Level level);

  const Config& config() const { return m_cfg; }

  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<SimVertexContainer> m_inputVertices{this, "InputVertices"};

  WriteDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_outputEvent{
      this, "OutputEvent"};
};

}  // namespace ActsExamples
