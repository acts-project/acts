// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepParticleWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <stdexcept>

#include "edm4hep/MCParticle.h"

namespace ActsExamples {

EDM4hepParticleWriter::EDM4hepParticleWriter(
    const EDM4hepParticleWriter::Config& config, Acts::Logging::Level lvl)
    : WriterT(config.inputParticles, "EDM4hepParticleWriter", lvl),
      m_cfg(config),
      m_writer(config.outputPath, &m_store) {
  ACTS_VERBOSE("Created output file " << config.outputPath);

  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }

  m_mcParticleCollection =
      &m_store.create<edm4hep::MCParticleCollection>(m_cfg.outputParticles);
  m_writer.registerForWrite(m_cfg.outputParticles);
}

ActsExamples::ProcessCode EDM4hepParticleWriter::endRun() {
  m_writer.finish();

  return ProcessCode::SUCCESS;
}

ProcessCode EDM4hepParticleWriter::writeT(
    const AlgorithmContext&, const SimParticleContainer& particles) {
  for (const auto& particle : particles) {
    auto p = m_mcParticleCollection->create();
    EDM4hepUtil::writeParticle(particle, p);
  }

  m_writer.writeEvent();
  m_store.clearCollections();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
