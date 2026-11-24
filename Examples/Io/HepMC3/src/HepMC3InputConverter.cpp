// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3InputConverter.hpp"

#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

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

HepMC3InputConverter::HepMC3InputConverter(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("HepMC3InputConverter", std::move(logger)), m_cfg(config) {
  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputVertices.initialize(m_cfg.outputVertices);
  m_inputEvent.initialize(m_cfg.inputEvent);
  if (m_cfg.outputDecayGraph.empty() != m_cfg.decayPdgCodes.empty()) {
    throw std::invalid_argument(
        "Specify both outputDecayGraph and decayPdgCodes");
  }
  if (!m_cfg.outputDecayGraph.empty()) {
    m_outputDecayGraph.initialize(m_cfg.outputDecayGraph);
  }
}

ProcessCode HepMC3InputConverter::execute(const AlgorithmContext& ctx) const {
  ACTS_DEBUG("HepMC3InputConverter::execute");
  const auto& genEvent = *m_inputEvent(ctx);
  ACTS_DEBUG("Have " << genEvent.vertices().size() << " vertices");
  ACTS_DEBUG("Have " << genEvent.particles().size() << " particles");
  ACTS_DEBUG("Have " << genEvent.event_number() << " event number");

  convertHepMC3ToInternalEdm(ctx, genEvent);

  return ProcessCode::SUCCESS;
}

namespace {

constexpr int kBeamParticleStatus = 4;
constexpr int kUndecayedParticleStatus = 1;

Acts::Vector4 convertPosition(const HepMC3::GenVertex& vertex) {
  const auto& vec = vertex.position();
  const double unit =
      vertex.parent_event()->length_unit() == HepMC3::Units::CM ? 1_cm : 1_mm;
  return unit * Acts::Vector4(vec.x(), vec.y(), vec.z(), vec.t());
}

Acts::Vector4 convertMomentum(const HepMC3::GenParticle& particle) {
  const auto& vec = particle.momentum();
  const double unit =
      particle.parent_event()->momentum_unit() == HepMC3::Units::MEV ? 1_MeV
                                                                     : 1_GeV;
  return unit * Acts::Vector4(vec.px(), vec.py(), vec.pz(), vec.e());
}

std::invalid_argument badDecay(const HepMC3::GenParticle& source,
                               const std::string& reason) {
  return std::invalid_argument("HepMC decay at particle " +
                               std::to_string(source.id()) + " (PDG " +
                               std::to_string(source.pid()) + "): " + reason);
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
    for (const auto& [idx, particleId] : Acts::enumerate(vertex.incoming)) {
      const auto& particle = findParticle(particleId);
      if (idx == 0) {
        ss << " I:     ";
      } else {
        ss << "        ";
      }
      ss << particle;
      ss << "\n";
    }

    for (const auto& [idx, particleId] : Acts::enumerate(vertex.outgoing)) {
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

ActsExamples::HeavyFlavourOrigin HepMC3InputConverter::deriveHeavyFlavourOrigin(
    const std::shared_ptr<const HepMC3::GenParticle>& particleToCheck) const {
  std::vector<std::shared_ptr<const HepMC3::GenParticle>> st;
  std::unordered_set<int> visited;
  st.push_back(particleToCheck);

  bool isFromCharm{false};
  while (!st.empty()) {
    const std::shared_ptr<const HepMC3::GenParticle> part = st.back();
    st.pop_back();

    if (!part) {
      continue;
    }

    if (visited.count(part->id()) != 0u) {
      continue;
    }
    visited.insert(part->id());

    int pdgCode = std::abs(part->pid());
    auto hadType = Acts::ParticleIdHelper::hadronType(
        static_cast<Acts::PdgParticle>(pdgCode));

    // --- beauty PDG IDs ---
    if (hadType == Acts::HadronType::BottomMeson ||
        hadType == Acts::HadronType::BottomBaryon ||
        hadType == Acts::HadronType::BBbarMeson ||
        (m_cfg.searchUpToHeavyFlavourQuark &&
         static_cast<ActsExamples::HeavyFlavourOrigin>(pdgCode) ==
             ActsExamples::HeavyFlavourOrigin::Bottom)) {
      return ActsExamples::HeavyFlavourOrigin::Bottom;
    }

    // --- charm PDG IDs ---
    if (hadType == Acts::HadronType::CharmedMeson ||
        hadType == Acts::HadronType::CharmedBaryon ||
        hadType == Acts::HadronType::CCbarMeson ||
        (m_cfg.searchUpToHeavyFlavourQuark &&
         static_cast<ActsExamples::HeavyFlavourOrigin>(pdgCode) ==
             ActsExamples::HeavyFlavourOrigin::Charm)) {
      // we do not return directly because
      // B -> D -> X should be tagged as from beauty and not charm
      isFromCharm = true;
    }

    // go to parents
    const auto& vtx = part->production_vertex();
    if (!vtx) {
      continue;
    }

    for (const auto& parent : vtx->particles_in()) {
      st.push_back(parent);
    }
  }

  if (isFromCharm) {
    return ActsExamples::HeavyFlavourOrigin::Charm;
  }

  return ActsExamples::HeavyFlavourOrigin::None;
}

std::optional<SimParticle> HepMC3InputConverter::makeParticle(
    const std::shared_ptr<const HepMC3::GenParticle>& source, SimBarcode id,
    const Acts::Vector4& position, bool prescribed) const {
  const Acts::PdgParticle pdg{source->pid()};
  auto massFromPdg = Acts::findMass(pdg);
  auto chargeFromPdg = Acts::findCharge(pdg);
  if (pdg == Acts::PdgParticle::eInvalid) {
    massFromPdg = chargeFromPdg = 0.;
  }
  if (prescribed) {
    // Geant4 transports neutrinos as massless. ACTS lists their masses as
    // unknown, so retain this convention only for the prescribed decay tree.
    const int absPdg = std::abs(source->pid());
    if (absPdg == 12 || absPdg == 14 || absPdg == 16) {
      massFromPdg = 0.;
    }
  }
  if (!massFromPdg || !chargeFromPdg) {
    if (prescribed) {
      throw badDecay(*source, "missing PDG mass or charge");
    }
    ACTS_ERROR("No mass or charge found for PDG ID " << pdg);
    throw std::bad_optional_access{};
  }
  const double mass = massFromPdg.value();
  const double charge = chargeFromPdg.value();
  const Acts::Vector4 mom = convertMomentum(*source);
  if (prescribed) {
    // Allow HepMC text rounding and small differences between PDG table
    // editions.
    if (pdg == Acts::PdgParticle::eInvalid || !mom.allFinite() ||
        !position.allFinite() || !std::isfinite(mass) ||
        !std::isfinite(charge) || mom[3] <= 0. ||
        std::abs(mom[3] - std::hypot(mass, mom.head<3>().norm())) >
            1e-6_GeV + 1e-6 * mom[3]) {
      throw badDecay(*source, "expected finite on-shell four-momentum");
    }
  } else if (std::isnan(mass) || std::isnan(charge)) {
    return std::nullopt;
  }

  SimParticle particle{id, pdg, charge, mass};
  particle.initialState().setPosition4(position);
  const double p = mom.head<3>().norm();
  if (p > 0.) {
    particle.initialState().setDirection(mom.head<3>() / p);
  }
  particle.initialState().setAbsoluteMomentum(p);
  particle.setOrigParticleIdx(source->id());
  particle.setHeavyFlavourOrigin(deriveHeavyFlavourOrigin(source));
  return particle;
}

void HepMC3InputConverter::addDecay(
    const HepMC3::GenParticle& source, const SimParticle& parent,
    SimParticleDecayGraph& graph, std::vector<bool>& seenVertices,
    SimBarcode::SubParticleId& nextSubparticle) const {
  const auto& end = *source.end_vertex();
  if (source.status() != 2 || end.particles_in().size() != 1 ||
      end.particles_out().size() < 2 ||
      seenVertices.at(std::abs(end.id()) - 1) || parent.mass() <= 0.) {
    throw badDecay(source,
                   "expected an unvisited single-parent decay with at least "
                   "two daughters");
  }
  seenVertices.at(std::abs(end.id()) - 1) = true;
  SimParticleDecay decay;
  decay.endPosition = convertPosition(end);
  const Acts::Vector4 flight = decay.endPosition - parent.fourPosition();
  const double ct = flight[3];
  const Acts::Vector3 expectedFlight = ct * parent.momentum() / parent.energy();
  if (!flight.allFinite() || ct < 0. ||
      (flight.head<3>() - expectedFlight).norm() > 1e-3_mm + 1e-5 * ct) {
    throw badDecay(source,
                   "inconsistent decay position, time and parent momentum");
  }
  decay.properTime = ct * parent.mass() / parent.energy();
  Acts::Vector4 daughterSum = Acts::Vector4::Zero();
  double daughterCharge = 0.;
  for (const auto& child : end.particles_out()) {
    if (parent.particleId().generation() ==
            std::numeric_limits<SimBarcode::GenerationId>::max() ||
        nextSubparticle ==
            std::numeric_limits<SimBarcode::SubParticleId>::max()) {
      throw badDecay(source, "decay tree exceeds the barcode range");
    }
    const SimBarcode id = parent.particleId().makeDescendant(++nextSubparticle);
    SimParticle daughter =
        makeParticle(child, id, decay.endPosition, true).value();
    daughter.setParentParticleId(parent.particleId());
    daughterSum += convertMomentum(*child);
    daughterCharge += daughter.charge();
    if (child->end_vertex() != nullptr) {
      if (std::ranges::find(m_cfg.decayPdgCodes, child->pid()) ==
          m_cfg.decayPdgCodes.end()) {
        throw badDecay(*child,
                       "nested decay requires its PDG in decayPdgCodes");
      }
      addDecay(*child, daughter, graph, seenVertices, nextSubparticle);
    } else if (child->status() != kUndecayedParticleStatus) {
      throw badDecay(*child, "terminal daughter must have status 1");
    }
    decay.daughters.push_back(std::move(daughter));
  }
  if (std::abs(daughterCharge - parent.charge()) > 1e-6) {
    throw badDecay(source, "daughter charges do not sum to the parent charge");
  }
  if ((daughterSum - convertMomentum(source)).cwiseAbs().maxCoeff() >
      1e-6_GeV + 1e-6 * parent.energy()) {
    throw badDecay(source,
                   "daughter four-momenta do not conserve four-momentum");
  }
  if (!graph.emplace(parent.particleId(), std::move(decay)).second) {
    throw badDecay(source, "duplicate generator parent barcode");
  }
}

void HepMC3InputConverter::handleVertex(
    const HepMC3::GenVertex& genVertex, SimVertex& vertex,
    std::vector<SimVertex>& vertices, std::vector<SimParticle>& particles,
    std::size_t& nSecondaryVertices, std::size_t& nParticles,
    std::vector<bool>& seenVertices, SimParticleDecayGraph* graph) const {
  for (const auto& particle : genVertex.particles_out()) {
    const bool prescribe =
        graph != nullptr && particle->status() == 2 &&
        std::ranges::find(m_cfg.decayPdgCodes, particle->pid()) !=
            m_cfg.decayPdgCodes.end();
    if (prescribe && !particle->end_vertex()) {
      throw badDecay(*particle, "status-2 parent is missing its decay vertex");
    }
    if (!prescribe && particle->end_vertex() != nullptr) {
      // This particle has an end vertex, we need to handle that vertex
      const auto& endVertex = *particle->end_vertex();
      if (seenVertices.at(std::abs(endVertex.id()) - 1)) {
        // We've seen this vertex before, we don't need to create it again
        continue;
      }
      seenVertices.at(std::abs(endVertex.id()) - 1) = true;

      const auto endVertexPosition = convertPosition(endVertex);

      ACTS_VERBOSE("Found secondary vertex at "
                   << endVertexPosition.transpose());
      double distance = (endVertexPosition.template head<3>() -
                         vertex.position4.template head<3>())
                            .squaredNorm();

      if (distance <=
              m_cfg.vertexSpatialThreshold * m_cfg.vertexSpatialThreshold &&
          m_cfg.mergeSecondaries) {
        handleVertex(endVertex, vertex, vertices, particles, nSecondaryVertices,
                     nParticles, seenVertices, graph);
      } else {
        // Over threshold, make a new vertex
        SimVertex secondaryVertex;
        nSecondaryVertices += 1;
        if (graph != nullptr &&
            nSecondaryVertices >
                std::numeric_limits<SimBarcode::SecondaryVertexId>::max()) {
          throw badDecay(*particle,
                         "secondary vertex exceeds the barcode range");
        }
        secondaryVertex.id =
            SimVertexBarcode{vertex.id}.withVertexSecondary(nSecondaryVertices);
        secondaryVertex.position4 = convertPosition(endVertex);

        handleVertex(endVertex, secondaryVertex, vertices, particles,
                     nSecondaryVertices, nParticles, seenVertices, graph);

        // Only keep the secondary vertex if it has outgoing particles
        if (!secondaryVertex.outgoing.empty()) {
          vertices.push_back(secondaryVertex);
        }
      }
    } else {
      if (!prescribe && particle->status() != kUndecayedParticleStatus) {
        ACTS_ERROR("Undecayed particle has status " << particle->status());
      }
      const SimBarcode particleId =
          SimBarcode()
              .withVertexPrimary(vertex.vertexId().vertexPrimary())
              .withVertexSecondary(vertex.vertexId().vertexSecondary())
              .withParticle(++nParticles);
      auto part = makeParticle(
          particle, particleId,
          prescribe ? convertPosition(genVertex) : vertex.position4, prescribe);
      if (!part) {
        continue;
      }
      if (prescribe) {
        SimBarcode::SubParticleId nextSubparticle = 0;
        addDecay(*particle, *part, *graph, seenVertices, nextSubparticle);
      }
      particles.push_back(std::move(*part));
      vertex.outgoing.insert(particleId);
    }
  }
}

void HepMC3InputConverter::convertHepMC3ToInternalEdm(
    const AlgorithmContext& ctx, const HepMC3::GenEvent& genEvent) const {
  ACTS_DEBUG("Converting HepMC3 event to internal EDM");

  ACTS_VERBOSE("Have " << genEvent.vertices().size() << " vertices");
  ACTS_VERBOSE("Have " << genEvent.particles().size() << " particles");

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

        auto position = convertPosition(*vertex);

        if (auto it = std::ranges::find_if(
                vertexClusters,
                [&](const auto& cluster) {
                  const auto clusterPosition = convertPosition(*cluster.at(0));
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

  SimParticleDecayGraph graph;
  std::vector<bool> seenVertices;
  seenVertices.resize(genEvent.vertices().size());

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
    ACTS_DEBUG("Reserving " << genEvent.vertices().size() / 2
                            << " vertices and " << nUndecayedParticles
                            << " particles");
    verticesUnordered.reserve(genEvent.vertices().size() / 2);
    particlesUnordered.reserve(nUndecayedParticles);

    for (auto& cluster : vertexClusters) {
      ACTS_VERBOSE("Primary vertex cluster at "
                   << convertPosition(*cluster.at(0)).transpose()
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
      if (!m_cfg.outputDecayGraph.empty() &&
          nPrimaryVertices >
              std::numeric_limits<SimBarcode::PrimaryVertexId>::max()) {
        throw std::overflow_error(
            "HepMC primary vertex exceeds the barcode range");
      }
      SimVertex primaryVertex;
      primaryVertex.id = SimVertexBarcode().withVertexPrimary(nPrimaryVertices);
      primaryVertex.position4 = convertPosition(*cluster.at(0));

      std::size_t nSecondaryVertices = 0;
      std::size_t nParticles = 0;

      for (auto& genVertex : cluster) {
        handleVertex(*genVertex, primaryVertex, verticesUnordered,
                     particlesUnordered, nSecondaryVertices, nParticles,
                     seenVertices,
                     m_cfg.outputDecayGraph.empty() ? nullptr : &graph);
      }
      verticesUnordered.push_back(primaryVertex);
    }
  }

  if (!m_cfg.outputDecayGraph.empty()) {
    const auto expected =
        std::ranges::count_if(genEvent.particles(), [&](const auto& source) {
          return source->status() == 2 &&
                 std::ranges::find(m_cfg.decayPdgCodes, source->pid()) !=
                     m_cfg.decayPdgCodes.end();
        });
    if (graph.size() != static_cast<std::size_t>(expected)) {
      throw std::invalid_argument(
          "HepMC contains selected decays unreachable from a primary vertex");
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
  if (!m_cfg.outputDecayGraph.empty()) {
    m_outputDecayGraph(ctx, std::move(graph));
  }
}

}  // namespace ActsExamples
