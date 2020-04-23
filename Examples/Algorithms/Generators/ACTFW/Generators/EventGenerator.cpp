// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Generators/EventGenerator.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>

#include "ACTFW/Framework/WhiteBoard.hpp"

FW::EventGenerator::EventGenerator(const Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EventGenerator", lvl)) {
  if (m_cfg.output.empty()) {
    throw std::invalid_argument("Missing output collection");
  }
  if (m_cfg.generators.empty()) {
    throw std::invalid_argument("No generators are configured");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers service");
  }
}

std::string FW::EventGenerator::name() const {
  return "EventGenerator";
}

std::pair<size_t, size_t> FW::EventGenerator::availableEvents() const {
  return {0u, SIZE_MAX};
}

FW::ProcessCode FW::EventGenerator::read(const AlgorithmContext& ctx) {
  std::vector<SimVertex> event;

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  // number of primary vertices within event
  size_t nPrimaryVertices = 0;
  // number of secondary vertices associated to a primary vertex
  size_t nSecondaryVertices = 0;
  // number of particles associated to a vertex (primary + secondary)
  size_t nParticlesVertex = 0;
  // total number of particles within event
  size_t nParticles = 0;

  for (size_t iGenerate = 0; iGenerate < m_cfg.generators.size(); ++iGenerate) {
    auto& generate = m_cfg.generators[iGenerate];

    // generate the number of primary vertices from this generator
    for (size_t n = generate.multiplicity(rng); 0 < n; --n) {
      nPrimaryVertices += 1;
      nSecondaryVertices = 0;
      nParticlesVertex = 0;

      // generate primary vertex position
      auto vertex = generate.vertex(rng);
      // generate associated process vertices
      // by convention the first process vertex should contain the
      // particles associated directly to the primary vertex itself.
      auto processVertices = generate.process(rng);

      for (auto& processVertex : processVertices) {
        nSecondaryVertices += 1;
        // shift the process vertex to the generated primary vertex position
        processVertex.position4 += vertex;

        auto updateParticleInPlace = [&](ActsFatras::Particle& particle) {
          // only set the primary vertex, leave everything else as-is
          // using the number of primary vertices as the index ensures
          // that barcode=0 is not used, since it is typically used elsewhere
          // to signify elements w/o an associated particle.
          const auto pid = ActsFatras::Barcode(particle.particleId())
                               .setVertexPrimary(nPrimaryVertices);
          // move particle to the vertex
          const auto pos4 = (vertex + particle.position4()).eval();
          // this changes the particle identity; must reassign.
          particle = particle.withParticleId(pid).setPosition4(pos4);
        };

        for (auto& particle : processVertex.incoming) {
          updateParticleInPlace(particle);
          nParticlesVertex += 1;
        }
        for (auto& particle : processVertex.outgoing) {
          updateParticleInPlace(particle);
          nParticlesVertex += 1;
        }
      }
      nParticles += nParticlesVertex;

      // append all process vertices to the full event
      std::move(processVertices.begin(), processVertices.end(),
                std::back_inserter(event));

      ACTS_VERBOSE("event=" << ctx.eventNumber << " generator=" << iGenerate
                            << " primary_vertex=" << nPrimaryVertices
                            << " n_secondary_vertices=" << nSecondaryVertices
                            << " n_particles=" << nParticlesVertex);
    }
  }
  // TODO should this reassign the vertex ids?
  // if not, what is the purpose? can it be removed altogether?
  if (m_cfg.shuffle) {
    std::shuffle(event.begin(), event.end(), rng);
  }

  ACTS_DEBUG("event=" << ctx.eventNumber
                      << " n_primary_vertices=" << nPrimaryVertices
                      << " n_secondary_vertices=" << event.size()
                      << " n_particles=" << nParticles);

  // move generated event to the store
  ctx.eventStore.add(m_cfg.output, std::move(event));
  return FW::ProcessCode::SUCCESS;
}
