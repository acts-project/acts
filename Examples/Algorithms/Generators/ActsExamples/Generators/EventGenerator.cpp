// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Generators/EventGenerator.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <ostream>
#include <span>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

// This is a hack to make HepMC3::Print::listing public
// It's pretty evil but should have no side-effects
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wkeyword-macro"
#endif
#define private public
#include <HepMC3/Print.h>
#undef private
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

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

namespace {
std::shared_ptr<HepMC3::GenEvent> mergeHepMC3Events(
    std::span<std::shared_ptr<HepMC3::GenEvent>> genEvents,
    const Acts::Logger& logger) {
  Acts::AveragingScopedTimer mergeTimer("Merging generator events", logger(),
                                        Acts::Logging::DEBUG);

  std::vector<std::shared_ptr<HepMC3::GenParticle>> particles;

  auto event = std::make_shared<HepMC3::GenEvent>();
  event->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

  for (auto& genEvent : genEvents) {
    auto sample = mergeTimer.sample();
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
      auto vertex = std::make_shared<HepMC3::GenVertex>(srcVertex->position());
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
  }

  return event;
}
}  // namespace

ProcessCode EventGenerator::read(const AlgorithmContext& ctx) {
  ACTS_VERBOSE("EventGenerator::read");
  std::vector<SimParticle> particlesUnordered;
  std::vector<SimVertex> verticesUnordered;

  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);

  std::vector<std::shared_ptr<HepMC3::GenEvent>> genEvents;

  std::size_t nPrimaryVertices = 0;
  ACTS_VERBOSE("Using " << m_cfg.generators.size() << " generators");
  {
    Acts::AveragingScopedTimer genTimer("Generating primary vertices", logger(),
                                        Acts::Logging::DEBUG);

    for (std::size_t iGenerate = 0; iGenerate < m_cfg.generators.size();
         ++iGenerate) {
      auto& generate = m_cfg.generators[iGenerate];

      // generate the primary vertices from this generator
      assert(generate.multiplicity != nullptr);
      std::size_t multiplicity = (*generate.multiplicity)(rng);
      ACTS_VERBOSE("Generating " << multiplicity << " primary vertices");
      for (std::size_t n = multiplicity; 0 < n; --n) {
        std::optional sample{genTimer.sample()};
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
        // Our internal time unit is ctau, so is HepMC3's, make sure we convert
        // to mm
        HepMC3::FourVector vtxPosHepMC(vertexPosition[Acts::eFreePos0] / 1_mm,
                                       vertexPosition[Acts::eFreePos1] / 1_mm,
                                       vertexPosition[Acts::eFreePos2] / 1_mm,
                                       vertexPosition[Acts::eFreeTime] / 1_mm);
        genEvent->shift_position_to(vtxPosHepMC);
        genEvent->add_attribute("acts",
                                std::make_shared<HepMC3::BoolAttribute>(true));

        if (m_cfg.printListing) {
          ACTS_VERBOSE("Generated event:\n"
                       << [&]() {
                            std::stringstream ss;
                            HepMC3::Print::listing(ss, *genEvent);
                            return ss.str();
                          }());
        }
        genEvents.push_back(genEvent);

        sample.reset();  // reset the gen timer
      }
    }
  }

  auto event = mergeHepMC3Events(genEvents, logger());
  event->set_event_number(ctx.eventNumber);

  ACTS_VERBOSE("Vertices size: " << event->vertices_size());
  if (m_cfg.printListing) {
    ACTS_VERBOSE("Final event:\n"
                 << [&]() {
                      std::stringstream ss;
                      HepMC3::Print::listing(ss, *event);
                      return ss.str();
                    }());
  }

  ACTS_DEBUG("event=" << ctx.eventNumber
                      << " n_primary_vertices=" << nPrimaryVertices
                      << " n_particles=" << event->particles_size());

  convertHepMC3ToInternalEdm(ctx, *event);

  if (m_outputEvent.isInitialized()) {
    m_outputEvent(ctx, std::move(event));
  }

  return ProcessCode::SUCCESS;
}

namespace {

constexpr int kBeamParticleStatus = 4;
constexpr int kUndecayedParticleStatus = 1;

Acts::Vector4 convertPosition(const HepMC3::FourVector& vec) {
  return Acts::Vector4(vec.x() * 1_mm, vec.y() * 1_mm, vec.z() * 1_mm,
                       vec.t() * 1_mm);
}

std::string printListing(const auto& vertices, const auto& particles) {
  auto findParticle = [&](SimBarcode particleId) {
    if (auto it = std::ranges::find_if(particles,
                                       [&](const auto& particle) {
                                         return particle.particleId() ==
                                                particleId;
                                       });
        it != particles.end()) {
      return *it;
    }
    throw std::runtime_error("Particle not found");
  };
  std::stringstream ss;

  for (const auto& vertex : vertices) {
    ss << "Vtx:    " << vertex.vertexId() << " at "
       << vertex.position4.transpose() << "\n";
    for (const auto& [idx, particleId] : enumerate(vertex.incoming)) {
      const auto& particle = findParticle(particleId);
      if (idx == 0) {
        ss << " I:     ";
      } else {
        ss << "        ";
      }
      ss << particle;
      ss << "\n";
    }

    for (const auto& [idx, particleId] : enumerate(vertex.outgoing)) {
      const auto& particle = findParticle(particleId);
      if (idx == 0) {
        ss << " O:     ";
      } else {
        ss << "        ";
      }
      ss << particle;
      ss << "\n";
    }
  }

  return ss.str();
};
}  // namespace

void EventGenerator::handleVertex(const HepMC3::GenVertex& genVertex,
                                  SimVertex& vertex,
                                  std::vector<SimVertex>& vertices,
                                  std::vector<SimParticle>& particles,
                                  std::size_t& nSecondaryVertices,
                                  std::size_t& nParticles,
                                  std::vector<bool>& seenVertices) {
  for (const auto& particle : genVertex.particles_out()) {
    if (particle->end_vertex() != nullptr) {
      // This particle has an end vertex, we need to handle that vertex
      const auto& endVertex = *particle->end_vertex();
      if (seenVertices.at(std::abs(endVertex.id()) - 1)) {
        // We've seen this vertex before, we don't need to create it again
        continue;
      }
      seenVertices.at(std::abs(endVertex.id()) - 1) = true;

      const auto endVertexPosition = convertPosition(endVertex.position());

      ACTS_VERBOSE("Found secondary vertex at "
                   << endVertexPosition.transpose());
      double distance = (endVertexPosition.template head<3>() -
                         vertex.position4.template head<3>())
                            .squaredNorm();

      if (distance <=
              m_cfg.vertexSpatialThreshold * m_cfg.vertexSpatialThreshold &&
          m_cfg.mergeSecondaries) {
        handleVertex(endVertex, vertex, vertices, particles, nSecondaryVertices,
                     nParticles, seenVertices);
      } else {
        // Over threshold, make a new vertex
        SimVertex secondaryVertex;
        nSecondaryVertices += 1;
        secondaryVertex.id =
            SimVertexBarcode{vertex.id}.setVertexSecondary(nSecondaryVertices);
        secondaryVertex.position4 = convertPosition(endVertex.position());

        handleVertex(endVertex, secondaryVertex, vertices, particles,
                     nSecondaryVertices, nParticles, seenVertices);

        // Only keep the secondary vertex if it has outgoing particles
        if (!secondaryVertex.outgoing.empty()) {
          vertices.push_back(secondaryVertex);
        }
      }
    } else {
      if (particle->status() != kUndecayedParticleStatus) {
        ACTS_ERROR("Undecayed particle has status "
                   << particle->status() << "(and not "
                   << kUndecayedParticleStatus << ")");
      }
      // This particle is a final state particle
      SimBarcode particleId{0u};
      nParticles += 1;
      particleId.setVertexPrimary(vertex.vertexId().vertexPrimary())
          .setVertexSecondary(vertex.vertexId().vertexSecondary())
          .setParticle(nParticles);

      Acts::PdgParticle pdg{particle->pdg_id()};
      double mass = Acts::findMass(pdg).value();
      double charge = Acts::findCharge(pdg).value();

      SimParticle simParticle{particleId, pdg, charge, mass};
      simParticle.initial().setPosition4(vertex.position4);

      const HepMC3::FourVector& genMomentum = particle->momentum();
      Acts::Vector3 momentum{genMomentum.px() * 1_GeV, genMomentum.py() * 1_GeV,
                             genMomentum.pz() * 1_GeV};
      double p = momentum.norm();

      simParticle.initial().setDirection(momentum.normalized());
      simParticle.initial().setAbsoluteMomentum(p);

      particles.push_back(simParticle);
      vertex.outgoing.insert(particleId);
    }
  }
}

void EventGenerator::convertHepMC3ToInternalEdm(
    const AlgorithmContext& ctx, const HepMC3::GenEvent& genEvent) {
  ACTS_DEBUG("Converting HepMC3 event to internal EDM");

  ACTS_VERBOSE("Have " << genEvent.vertices_size() << " vertices");
  ACTS_VERBOSE("Have " << genEvent.particles_size() << " particles");

  using VertexCluster = std::vector<HepMC3::ConstGenVertexPtr>;
  std::vector<VertexCluster> vertexClusters;

  ACTS_VERBOSE("Finding primary vertex clusters with threshold "
               << m_cfg.primaryVertexSpatialThreshold);

  // Find all vertices whose incoming particles are either all beam particles or
  // that don't have any incoming particles
  {
    Acts::ScopedTimer timer("Finding primary vertex clusters", logger(),
                            Acts::Logging::DEBUG);
    // @TODO: Maybe to it the other way round: go through particles and check if we already have that vertex?
    for (const auto& vertex : genEvent.vertices()) {
      if (vertex->particles_in().empty() ||
          std::ranges::all_of(vertex->particles_in(), [](const auto& particle) {
            return particle->status() == kBeamParticleStatus;
          })) {
        // Check if primary vertex is within tolerance of an existing primary
        // vertex

        ACTS_VERBOSE("Found primary vertex candidate at:\n"
                     << [&]() {
                          std::stringstream ss;
                          HepMC3::Print::listing(ss, vertex);
                          return ss.str();
                        }());

        auto position = convertPosition(vertex->position());

        if (auto it = std::ranges::find_if(
                vertexClusters,
                [&](const auto& cluster) {
                  const auto clusterPosition =
                      convertPosition(cluster.at(0)->position());
                  return (position - clusterPosition)
                             .template head<3>()
                             .cwiseAbs()
                             .maxCoeff() < m_cfg.primaryVertexSpatialThreshold;
                });
            it != vertexClusters.end() && m_cfg.mergePrimaries) {
          // Add the vertex to the cluster
          it->push_back(vertex);
        } else {
          // Make a new cluster
          vertexClusters.push_back({vertex});
        }
      }
    }
  }

  ACTS_DEBUG("Found " << vertexClusters.size()
                      << " primary vertex clusters (~ primary vertices)");

  std::vector<SimVertex> verticesUnordered;
  std::vector<SimParticle> particlesUnordered;

  std::vector<bool> seenVertices;
  seenVertices.resize(genEvent.vertices_size());

  for (const auto& cluster : vertexClusters) {
    for (const auto& vertex : cluster) {
      seenVertices.at(std::abs(vertex->id()) - 1) = true;
    }
  }

  std::size_t nPrimaryVertices = 0;
  {
    Acts::ScopedTimer timer("Converting HepMC3 vertices to internal EDM",
                            logger(), Acts::Logging::DEBUG);

    std::size_t nUndecayedParticles = 0;
    for (const auto& particle : genEvent.particles()) {
      // if (particle->status() == kUndecayedParticleStatus) {
      if (particle->end_vertex() == nullptr) {
        nUndecayedParticles += 1;
      }
    }
    ACTS_DEBUG("Found " << nUndecayedParticles
                        << " undecayed particles in the event");
    ACTS_DEBUG("Reserving " << genEvent.vertices_size() / 2 << " vertices and "
                            << nUndecayedParticles << " particles");
    verticesUnordered.reserve(genEvent.vertices_size() / 2);
    particlesUnordered.reserve(nUndecayedParticles);

    for (auto& cluster : vertexClusters) {
      ACTS_VERBOSE("Primary vertex cluster at "
                   << convertPosition(cluster.at(0)->position()).transpose()
                   << " containing " << cluster.size() << " vertices:\n"
                   <<
                   [&]() {
                     std::stringstream ss;
                     for (auto& vertex : cluster) {
                       HepMC3::Print::listing(ss, vertex);
                     }
                     return ss.str();
                   }()
                   << "--------------");

      nPrimaryVertices += 1;
      SimVertex primaryVertex;
      primaryVertex.id = SimVertexBarcode{}.setVertexPrimary(nPrimaryVertices);
      primaryVertex.position4 = convertPosition(cluster.at(0)->position());

      std::size_t nSecondaryVertices = 0;
      std::size_t nParticles = 0;

      for (auto& genVertex : cluster) {
        handleVertex(*genVertex, primaryVertex, verticesUnordered,
                     particlesUnordered, nSecondaryVertices, nParticles,
                     seenVertices);
      }
      verticesUnordered.push_back(primaryVertex);
    }
  }

  ACTS_DEBUG("Converted " << particlesUnordered.size() << " particles and "
                          << verticesUnordered.size() << " vertices");
  // for(std::size_t i =0;i<)

  if (m_cfg.printListing) {
    ACTS_VERBOSE("Converted event record:\n"
                 << printListing(verticesUnordered, particlesUnordered));
  }

  {
    Acts::ScopedTimer timer("Sorting particles and vertices", logger(),
                            Acts::Logging::DEBUG);

    std::ranges::sort(particlesUnordered, [](const auto& a, const auto& b) {
      return a.particleId() < b.particleId();
    });
    std::ranges::sort(verticesUnordered, [](const auto& a, const auto& b) {
      return a.vertexId() < b.vertexId();
    });
  }

  if (m_cfg.checkConsistency) {
    ACTS_DEBUG("Checking consistency of particles");
    auto equalParticleIds = [](const auto& a, const auto& b) {
      return a.particleId() == b.particleId();
    };
    if (auto it =
            std::ranges::adjacent_find(particlesUnordered, equalParticleIds);
        it != particlesUnordered.end()) {
      ACTS_ERROR("Found duplicate particle id " << it->particleId());
    }
    auto equalVertexIds = [](const auto& a, const auto& b) {
      return a.vertexId() == b.vertexId();
    };
    if (auto it = std::ranges::adjacent_find(verticesUnordered, equalVertexIds);
        it != verticesUnordered.end()) {
      ACTS_ERROR("Found duplicate vertex id " << it->vertexId());
    }
  }

  SimParticleContainer particles{particlesUnordered.begin(),
                                 particlesUnordered.end()};
  SimVertexContainer vertices{verticesUnordered.begin(),
                              verticesUnordered.end()};

  // std::vector<SimVertexBarcode> vertexBarcodes;
  // vertexBarcodes.resize(genEvent.vertices_size());

  // for(const auto& inParticle : genEvent.particles()) {
  //   // According to https://arxiv.org/pdf/1912.08005 status 1 is "Undecayed
  //   physical particle" if (inParticle->status() != 1) {
  //     continue;
  //   }

  //   SimBarcode particleId{0u};
  //   particleId.setParticle(particles.size() + 1);

  //   const auto* productionVertex = inParticle->production_vertex().get();
  //   if (productionVertex != nullptr) {
  //     const int vertexIndex = productionVertex->id();
  //     SimVertexBarcode vertexBarcode = vertexBarcodes.at(vertexIndex);
  //     if (vertexBarcode == SimVertexBarcode{0u}) {
  //       // We have not seen this vertex before, create a new one
  //       SimVertex vertex;
  //     }
  //   }
  // }

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

  // move generated event to the store
  m_outputParticles(ctx, std::move(particles));
  m_outputVertices(ctx, std::move(vertices));
}

}  // namespace ActsExamples
