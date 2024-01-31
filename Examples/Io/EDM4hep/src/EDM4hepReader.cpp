// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples//EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <algorithm>
#include <iomanip>
#include <stdexcept>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/Frame.h>
#include <podio/ObjectID.h>

namespace ActsExamples {

EDM4hepReader::EDM4hepReader(const Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepParticleReader", level)) {
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  m_reader.openFile(m_cfg.inputPath);

  m_eventsRange = std::make_pair(0, m_reader.getEntries("events"));

  m_outputParticles.initialize(m_cfg.outputParticles);
}

std::string EDM4hepReader::name() const {
  return "EDM4hepReader";
}

std::pair<std::size_t, std::size_t> EDM4hepReader::availableEvents() const {
  return m_eventsRange;
}

namespace {
std::string vid(unsigned int vtx) {
  return "V" + std::to_string(vtx);
}

std::string pid(const SimParticle& particle) {
  return "P" + std::to_string(particle.particleId().value());
};

std::string plabel(const SimParticle& particle) {
  using namespace Acts::UnitLiterals;
  std::stringstream ss;
  ss << particle.pdg() << "\\n(" << particle.particleId() << ")\\n"
     << "p=" << std::setprecision(3) << particle.absoluteMomentum() / 1_GeV
     << " GeV";
  return ss.str();
};

}  // namespace

void EDM4hepReader::graphviz(
    std::ostream& os, const SimParticleContainer::sequence_type& particles,
    const ParentRelationship& parents) const {
  os << "digraph Event {\n";

  std::set<unsigned int> primaryVertices;

  for (const auto& particle : particles) {
    if (particle.particleId().generation() == 0) {
      primaryVertices.insert(particle.particleId().vertexPrimary());

      os << vid(particle.particleId().vertexPrimary()) << " -> "
         << pid(particle) << ";\n";
    }

    os << pid(particle) << " [label=\"" << plabel(particle) << "\"];\n";
  }

  for (const auto [childIdx, parentIdx] : parents) {
    const auto& child = particles[childIdx];
    const auto& parent = particles[parentIdx];
    os << pid(parent) << " -> " << pid(child);

    if (parent.particleId().vertexSecondary() ==
        child.particleId().vertexSecondary()) {
      os << " [style=dashed]";
    }

    os << ";\n";
  }

  for (unsigned int vtx : primaryVertices) {
    os << vid(vtx) << " [label=\"PV" << vtx << "\"];\n";
  }

  os << "}";
}

ProcessCode EDM4hepReader::read(const AlgorithmContext& ctx) {
  podio::Frame frame = m_reader.readEntry("events", ctx.eventNumber);
  const auto& mcParticleCollection =
      frame.get<edm4hep::MCParticleCollection>(m_cfg.inputParticles);

  ACTS_DEBUG("Reading EDM4hep inputs");

  SimParticleContainer::sequence_type unordered;

  // Read particles from the input file
  // Find particles without parents and group them by vtx position to find
  // primary vertices
  std::vector<std::pair<Acts::Vector3, std::vector<edm4hep::MCParticle>>>
      primaryVertices;
  for (const auto& particle : mcParticleCollection) {
    if (particle.parents_size() > 0) {
      // not a primary vertex
      continue;
    }
    const auto& vtx = particle.getVertex();
    Acts::Vector3 vtxPos = {vtx[0], vtx[1], vtx[2]};
    vtxPos /= Acts::UnitConstants::mm;

    // linear search for vector
    auto it = std::find_if(
        primaryVertices.begin(), primaryVertices.end(),
        [&vtxPos](
            const std::pair<Acts::Vector3, std::vector<edm4hep::MCParticle>>&
                pair) { return pair.first == vtxPos; });

    if (it == primaryVertices.end()) {
      ACTS_DEBUG("Found primary vertex at " << vtx.x << ", " << vtx.y << ", "
                                            << vtx.z);
      primaryVertices.push_back({vtxPos, {particle}});
    } else {
      it->second.push_back(particle);
    }
  }

  ACTS_DEBUG("Found " << primaryVertices.size() << " primary vertices");

  // key: child, value: parent
  ParentRelationship parentRelationship;

  // key: input particle address, value: index in the unordered particle
  // container
  std::unordered_map<int, std::size_t> edm4hepParticleMap;

  std::size_t nPrimaryVertices = 0;
  // Walk the particle tree
  for (const auto& [vtxPos, particles] : primaryVertices) {
    nPrimaryVertices += 1;
    ACTS_DEBUG("Walking particle tree for primary vertex at "
               << vtxPos.x() << ", " << vtxPos.y() << ", " << vtxPos.z());
    std::size_t nParticles = 0;
    std::size_t nSecondaryVertices = 0;
    std::size_t maxGen = 0;
    auto startSize = unordered.size();
    for (const auto& inParticle : particles) {
      nParticles += 1;
      SimParticle particle{EDM4hepUtil::readParticle(inParticle)};
      particle.setParticleId(SimBarcode{}
                                 .setParticle(nParticles)
                                 .setVertexPrimary(nPrimaryVertices));
      ACTS_VERBOSE("+ add particle " << particle);
      ACTS_VERBOSE("  - at " << particle.position().transpose());
      ACTS_VERBOSE("  - createdInSim: " << inParticle.isCreatedInSimulation());
      ACTS_VERBOSE("  - vertexIsNotEndpointOfParent: "
                   << inParticle.vertexIsNotEndpointOfParent());
      ACTS_VERBOSE("  - isStopped: " << inParticle.isStopped());
      ACTS_VERBOSE("  - endpoint: " << inParticle.getEndpoint().x << ", "
                                    << inParticle.getEndpoint().y << ", "
                                    << inParticle.getEndpoint().z);
      const auto pid = particle.particleId();
      unordered.push_back(std::move(particle));
      edm4hepParticleMap[inParticle.getObjectID().index] = unordered.size() - 1;
      processChildren(inParticle, pid, unordered, parentRelationship,
                      edm4hepParticleMap, nSecondaryVertices, maxGen);
    }
    ACTS_VERBOSE("Primary vertex complete, produced "
                 << (unordered.size() - startSize) << " particles and "
                 << nSecondaryVertices << " secondary vertices in " << maxGen
                 << " generations");
    setSubParticleIds(std::next(unordered.begin(), startSize), unordered.end());
  }

  ACTS_DEBUG("Found " << unordered.size() << " particles");

  // Write ordered particles container to the EventStore
  SimParticleContainer particles;
  particles.insert(unordered.begin(), unordered.end());

  if (!m_cfg.graphvizOutput.empty()) {
    std::string path = perEventFilepath(m_cfg.graphvizOutput, "particles.dot",
                                        ctx.eventNumber);
    std::ofstream dot(path);
    graphviz(dot, unordered, parentRelationship);
  }

  ACTS_DEBUG("Reading sim hits from " << m_cfg.inputSimHits.size()
                                      << " sim hit collections");
  for (const auto& name : m_cfg.inputSimHits) {
    const auto& inputHits = frame.get<edm4hep::SimTrackerHitCollection>(name);

    for (const auto& hit : inputHits) {
      auto simHit = EDM4hepUtil::readSimHit(hit, [&](const auto& inParticle) {
        ACTS_VERBOSE("SimHit has source particle: "
                     << hit.getMCParticle().getObjectID().index);
        auto it = edm4hepParticleMap.find(inParticle.getObjectID().index);
        if (it == edm4hepParticleMap.end()) {
          ACTS_ERROR("SimHit has source particle that we did not see before");
          return SimBarcode{};
        }
        const auto& particle = unordered.at(it->second);
        ACTS_VERBOSE("- " << inParticle.getObjectID().index << " -> "
                          << particle.particleId());
        return particle.particleId();
      });
    }
  }

  m_outputParticles(ctx, std::move(particles));

  return ProcessCode::SUCCESS;
}

void EDM4hepReader::processChildren(
    const edm4hep::MCParticle& inParticle, SimBarcode parentId,
    SimParticleContainer::sequence_type& particles,
    ParentRelationship& parentRelationship,
    std::unordered_map<int, std::size_t>& particleMap,
    std::size_t& nSecondaryVertices, std::size_t& maxGen) const {
  constexpr auto indent = [&](std::size_t n) {
    std::string result;
    for (std::size_t i = 0; i < n; ++i) {
      result += "  ";
    }
    return result;
  };

  const std::size_t gen = parentId.generation();
  maxGen = std::max(maxGen, gen);

  ACTS_VERBOSE(indent(gen) << "  - processing daughters for input particle "
                           << inParticle.id());
  ACTS_VERBOSE(indent(gen) << "    -> found " << inParticle.daughters_size()
                           << " daughter(s)");

  bool parentDecayed =
      std::any_of(inParticle.daughters_begin(), inParticle.daughters_end(),
                  [](const edm4hep::MCParticle& daughter) {
                    return !daughter.vertexIsNotEndpointOfParent();
                  });
  std::size_t secondaryVertex = 0;
  if (parentDecayed) {
    ACTS_VERBOSE(indent(gen) << "    -> parent decays");
    secondaryVertex = ++nSecondaryVertices;
  }

  std::size_t parentIndex = particles.size() - 1;

  std::size_t nParticles = 0;
  for (const auto& daughter : inParticle.getDaughters()) {
    SimParticle particle = EDM4hepUtil::readParticle(daughter);

    auto pid = parentId.makeDescendant(nParticles++);
    if (daughter.vertexIsNotEndpointOfParent()) {
      // incoming particle survived, interaction via descendant
    } else {
      // incoming particle decayed
      pid = pid.setVertexSecondary(secondaryVertex);
    }
    particle.setParticleId(pid);

    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "+ add particle " << particle);
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "  - generation: " << particle.particleId().generation());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "  - at " << particle.position().transpose());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - createdInSim: "
                 << daughter.isCreatedInSimulation());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - vertexIsNotEndpointOfParent: "
                 << daughter.vertexIsNotEndpointOfParent());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - isStopped: " << daughter.isStopped());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - endpoint: " << daughter.getEndpoint().x << ", "
                 << daughter.getEndpoint().y << ", "
                 << daughter.getEndpoint().z);

    particles.push_back(std::move(particle));
    particleMap[daughter.getObjectID().index] = particles.size() - 1;
    parentRelationship[particles.size() - 1] = parentIndex;
    processChildren(daughter, pid, particles, parentRelationship, particleMap,
                    nSecondaryVertices, maxGen);
  }
}

void EDM4hepReader::setSubParticleIds(
    SimParticleContainer::sequence_type::iterator begin,
    SimParticleContainer::sequence_type::iterator end) {
  std::vector<std::size_t> numByGeneration;
  numByGeneration.reserve(10);

  for (auto it = begin; it != end; ++it) {
    auto& particle = *it;
    const auto pid = particle.particleId();
    if (pid.generation() >= numByGeneration.size()) {
      numByGeneration.resize(pid.generation() + 1, 0);
    }
    unsigned int nextSubParticle = numByGeneration[pid.generation()]++;

    auto newPid = particle.particleId().setSubParticle(nextSubParticle);
    particle.setParticleId(newPid);
  }
}

}  // namespace ActsExamples
