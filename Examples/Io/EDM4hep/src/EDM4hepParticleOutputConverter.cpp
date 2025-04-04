// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepParticleOutputConverter.hpp"

#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"

#include <stdexcept>

#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepParticleOutputConverter::EDM4hepParticleOutputConverter(
    const EDM4hepParticleOutputConverter::Config& cfg, Acts::Logging::Level lvl)
    : EDM4hepOutputConverter("EDM4hepParticleOutputConverter", lvl),
      m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }

  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing particles output collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputParticles.initialize(m_cfg.outputParticles);
}

ProcessCode EDM4hepParticleOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const SimParticleContainer particles = m_inputParticles(ctx);

  edm4hep::MCParticleCollection mcParticleCollection;

  for (const auto& particle : particles) {
    auto p = mcParticleCollection->create();
    EDM4hepUtil::writeParticle(particle, p);
  }

  m_outputParticles(ctx, std::move(mcParticleCollection));

  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepParticleOutputConverter::collections() const {
  return {m_cfg.outputParticles};
}

}  // namespace ActsExamples
