// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "PrintParticles.hpp"

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"

FW::PrintParticles::PrintParticles(const Config& cfg, Acts::Logging::Level lvl)
    : BareAlgorithm("PrintParticles", lvl), m_cfg(cfg) {}

FW::ProcessCode FW::PrintParticles::execute(
    const FW::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  for (const auto& particle : particles) {
    ACTS_INFO(
        "evt=" << ctx.eventNumber << " pid=" << particle.particleId()
               << " pdg=" << particle.pdg() << " phi="
               << Acts::VectorHelpers::phi(particle.unitDirection()) / 1_degree
               << "deg"
               << " eta=" << Acts::VectorHelpers::eta(particle.unitDirection())
               << " p=" << particle.absMomentum() / 1_GeV << "GeV");
  }
  return ProcessCode::SUCCESS;
}
