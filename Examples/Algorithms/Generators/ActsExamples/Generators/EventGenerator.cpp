// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/EventGenerator.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Print.h>

using namespace Acts::UnitLiterals;

namespace ActsExamples {

EventGenerator::EventGenerator(const Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg), m_logger(Acts::getDefaultLogger("EventGenerator", lvl)) {
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }
  if (m_cfg.outputVertices.empty()) {
    throw std::invalid_argument("Missing output vertices collection");
  }
  if (m_cfg.generators.empty()) {
    throw std::invalid_argument("No generators are configured");
  }

  for (auto& generator : m_cfg.generators) {
    if (generator.multiplicity == nullptr) {
      throw std::invalid_argument("Missing multiplicity generator");
    }
    if (generator.vertex == nullptr) {
      throw std::invalid_argument("Missing vertex generator");
    }
    if (generator.particles == nullptr) {
      throw std::invalid_argument("Missing particles generator");
    }
  }

  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers service");
  }

  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputVertices.initialize(m_cfg.outputVertices);
  m_outputEvent.maybeInitialize(m_cfg.outputEvent);
}

std::string EventGenerator::name() const {
  return "EventGenerator";
}

std::pair<std::size_t, std::size_t> EventGenerator::availableEvents() const {
  return {0u, std::numeric_limits<std::size_t>::max()};
}

ProcessCode EventGenerator::read(const AlgorithmContext& ctx) {
  ACTS_VERBOSE("EventGenerator::read");
  std::vector<SimParticle> particlesUnordered;
  std::vector<SimVertex> verticesUnordered;

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  auto event = std::make_shared<HepMC3::GenEvent>();
  event->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
  event->set_event_number(ctx.eventNumber);

  std::size_t nPrimaryVertices = 0;
  ACTS_VERBOSE("Using " << m_cfg.generators.size() << " generators");
  for (std::size_t iGenerate = 0; iGenerate < m_cfg.generators.size();
       ++iGenerate) {
    auto& generate = m_cfg.generators[iGenerate];

    std::vector<std::shared_ptr<HepMC3::GenParticle>> particles;

    // generate the primary vertices from this generator
    assert(generate.multiplicity != nullptr);
    std::size_t multiplicity = (*generate.multiplicity)(rng);
    ACTS_VERBOSE("Generating " << multiplicity << " primary vertices");
    for (std::size_t n = multiplicity; 0 < n; --n) {
      nPrimaryVertices += 1;

      // generate primary vertex position
      assert(generate.vertex != nullptr);
      auto vertexPosition = (*generate.vertex)(rng);
      ACTS_VERBOSE("Generate vertex at " << vertexPosition.transpose());

      // generate particles associated to this vertex
      assert(generate.particles != nullptr);
      auto genEvent = (*generate.particles)(rng);
      if (genEvent->length_unit() != HepMC3::Units::MM) {
        throw std::runtime_error("EventGenerator: length unit is not mm");
      }
      if (genEvent->momentum_unit() != HepMC3::Units::GEV) {
        throw std::runtime_error("EventGenerator: momentum unit is not GeV");
      }

      ACTS_VERBOSE("Shifting event to " << vertexPosition.transpose());
      // Our internal time unit is ctau, so is HepMC3's, make sure we convert to
      // mm
      HepMC3::FourVector vtxPosHepMC(vertexPosition[Acts::eFreePos0] / 1_mm,
                                     vertexPosition[Acts::eFreePos1] / 1_mm,
                                     vertexPosition[Acts::eFreePos2] / 1_mm,
                                     vertexPosition[Acts::eFreeTime] / 1_mm);
      genEvent->shift_position_to(vtxPosHepMC);
      genEvent->add_attribute("acts",
                              std::make_shared<HepMC3::BoolAttribute>(true));

      ACTS_VERBOSE("Generated event:\n"
                   << [&]() {
                        std::stringstream ss;
                        HepMC3::Print::content(ss, *genEvent);
                        return ss.str();
                      }());

      particles.clear();
      particles.reserve(genEvent->particles_size());

      auto copyAttributes = [&](const auto& src, auto& dst) {
        for (auto& attr : src.attribute_names()) {
          auto value = src.attribute_as_string(attr);
          dst.add_attribute(attr,
                            std::make_shared<HepMC3::StringAttribute>(value));
        }
      };

      copyAttributes(*genEvent, *event);

      // Add to combined event
      for (auto& srcParticle : genEvent->particles()) {
        if (srcParticle->id() - 1 != static_cast<int>(particles.size())) {
          throw std::runtime_error("Particle id is not consecutive");
        }
        auto particle = std::make_shared<HepMC3::GenParticle>();
        particle->set_momentum(srcParticle->momentum());
        particle->set_generated_mass(srcParticle->generated_mass());
        particle->set_pid(srcParticle->pid());
        particle->set_status(srcParticle->status());

        particles.push_back(particle);
        event->add_particle(particle);

        copyAttributes(*srcParticle, *particle);
      }

      for (auto& srcVertex : genEvent->vertices()) {
        auto vertex =
            std::make_shared<HepMC3::GenVertex>(srcVertex->position());
        vertex->set_status(srcVertex->status());

        event->add_vertex(vertex);

        copyAttributes(*srcVertex, *vertex);

        for (auto& srcParticle : srcVertex->particles_in()) {
          auto& particle = particles.at(srcParticle->id() - 1);
          vertex->add_particle_in(particle);
        }
        for (auto& srcParticle : srcVertex->particles_out()) {
          auto& particle = particles.at(srcParticle->id() - 1);
          vertex->add_particle_out(particle);
        }
      }

      // auto updateParticleInPlace = [&](SimParticle& particle) {
      //   // only set the primary vertex, leave everything else as-is
      //   // using the number of primary vertices as the index ensures
      //   // that barcode=0 is not used, since it is used elsewhere
      //   // to signify elements w/o an associated particle.
      //   const auto pid = SimBarcode{particle.particleId()}.setVertexPrimary(
      //       nPrimaryVertices);
      //   // move particle to the vertex
      //   const auto pos4 = (vertexPosition + particle.fourPosition()).eval();
      //   ACTS_VERBOSE(" - particle at " << pos4.transpose());
      //   // `withParticleId` returns a copy because it changes the identity
      //   particle = particle.withParticleId(pid);
      //   particle.initial().setPosition4(pos4);
      //   particle.final().setPosition4(pos4);
      // };
      // for (auto& vertexParticle : newParticles) {
      //   updateParticleInPlace(vertexParticle);
      // }

      // auto updateVertexInPlace = [&](SimVertex& vertex) {
      //   // only set the primary vertex, leave everything else as-is
      //   // using the number of primary vertices as the index ensures
      //   // that barcode=0 is not used, since it is used elsewhere
      //   // to signify elements w/o an associated particle.
      //   vertex.id = SimVertexBarcode{vertex.vertexId()}.setVertexPrimary(
      //       nPrimaryVertices);
      //   // move vertex
      //   const auto pos4 = (vertexPosition + vertex.position4).eval();
      //   ACTS_VERBOSE(" - vertex at " << pos4.transpose());
      //   vertex.position4 = pos4;

      //   for (auto& particleId : vertex.incoming) {
      //     particleId =
      //         SimBarcode{particleId}.setVertexPrimary(nPrimaryVertices);
      //   }
      //   for (auto& particleId : vertex.outgoing) {
      //     particleId =
      //         SimBarcode{particleId}.setVertexPrimary(nPrimaryVertices);
      //   }
      // };
      // for (auto& vertex : newVertices) {
      //   updateVertexInPlace(vertex);
      // }

      // ACTS_VERBOSE("event=" << ctx.eventNumber << " generator=" << iGenerate
      //                       << " primary_vertex=" << nPrimaryVertices
      //                       << " n_particles=" << newParticles.size());

      // particles.merge(std::move(newParticles));
      // vertices.merge(std::move(newVertices));
    }
  }

  ACTS_VERBOSE("Vertices size: " << event->vertices_size());
  ACTS_VERBOSE("Final event:\n"
               << [&]() {
                    std::stringstream ss;
                    HepMC3::Print::content(ss, *event);
                    return ss.str();
                  }());

  SimParticleContainer particles;
  SimVertexContainer vertices;

  ACTS_DEBUG("event=" << ctx.eventNumber
                      << " n_primary_vertices=" << nPrimaryVertices
                      << " n_particles=" << particles.size());

  // move generated event to the store
  m_outputParticles(ctx, std::move(particles));
  m_outputVertices(ctx, std::move(vertices));

  if (m_outputEvent.isInitialized()) {
    m_outputEvent(ctx, std::move(event));
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
