// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3InputConverter.hpp"

#include "Acts/Utilities/ScopedTimer.hpp"

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

HepMC3InputConverter::HepMC3InputConverter(const Config& config,
                                           Acts::Logging::Level level)
    : IAlgorithm("HepMC3InputConverter", level), m_cfg(config) {
  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputVertices.initialize(m_cfg.outputVertices);
  m_inputEvent.initialize(m_cfg.inputEvent);
}

ProcessCode HepMC3InputConverter::execute(const AlgorithmContext& ctx) const {
  ACTS_DEBUG("HepMC3InputConverter::execute");
  const auto& genEvent = *m_inputEvent(ctx);
  ACTS_DEBUG("Have " << genEvent.vertices_size() << " vertices");
  ACTS_DEBUG("Have " << genEvent.particles_size() << " particles");
  ACTS_DEBUG("Have " << genEvent.event_number() << " event number");

  convertHepMC3ToInternalEdm(ctx, genEvent);

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

void HepMC3InputConverter::handleVertex(const HepMC3::GenVertex& genVertex,
                                        SimVertex& vertex,
                                        std::vector<SimVertex>& vertices,
                                        std::vector<SimParticle>& particles,
                                        std::size_t& nSecondaryVertices,
                                        std::size_t& nParticles,
                                        std::vector<bool>& seenVertices) const {
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
      double mass = 0.0;
      double charge = 0.0;

      if (pdg != Acts::PdgParticle::eInvalid) {
        if (auto massOpt = Acts::findMass(pdg); massOpt.has_value()) {
          mass = massOpt.value();
        } else {
          ACTS_ERROR("No mass found for PDG ID " << pdg);
          throw std::bad_optional_access{};
        }

        if (auto chargeOpt = Acts::findCharge(pdg); chargeOpt.has_value()) {
          charge = chargeOpt.value();
        } else {
          ACTS_ERROR("No charge found for PDG ID " << pdg);
          throw std::bad_optional_access{};
        }
      }

      if (std::isnan(mass) || std::isnan(charge)) {
        // This particle does not have a mass or a charge based on our data
        // table:
        // => We can't handle it
        continue;
      }

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

void HepMC3InputConverter::convertHepMC3ToInternalEdm(
    const AlgorithmContext& ctx, const HepMC3::GenEvent& genEvent) const {
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

  // move generated event to the store
  m_outputParticles(ctx, std::move(particles));
  m_outputVertices(ctx, std::move(vertices));
}

}  // namespace ActsExamples
