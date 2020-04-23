// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/FlattenEvent.hpp"

#include <vector>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

FW::FlattenEvent::FlattenEvent(const Config& cfg, Acts::Logging::Level lvl)
    : BareAlgorithm("FlattenEvent", lvl), m_cfg(cfg) {
  if (m_cfg.inputEvent.empty()) {
    throw std::invalid_argument("Missing input event collection");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }
}

FW::ProcessCode FW::FlattenEvent::execute(const AlgorithmContext& ctx) const {
  // setup input and output containers
  const auto& event =
      ctx.eventStore.get<std::vector<SimVertex>>(m_cfg.inputEvent);
  SimParticleContainer::sequence_type unsortedParticles;

  // extract particles
  for (const auto& vertex : event) {
    for (const auto& particle : vertex.outgoing) {
      unsortedParticles.emplace_back(particle);
    }
  }
  unsortedParticles.shrink_to_fit();

  // re-establish ordering by barcode
  SimParticleContainer particles;
  particles.adopt_sequence(std::move(unsortedParticles));

  ctx.eventStore.add(m_cfg.outputParticles, std::move(particles));
  return ProcessCode::SUCCESS;
}
