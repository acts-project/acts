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
#include <stdexcept>

#include <edm4hep/MCParticle.h>
#include <podio/Frame.h>

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
      processChildren(inParticle, pid, unordered, parentRelationship,
                      nSecondaryVertices, maxGen);
    }
    ACTS_VERBOSE("Primary vertex complete, produced "
                 << (unordered.size() - startSize) << " particles and "
                 << nSecondaryVertices << " secondary vertices in " << maxGen
                 << " generations");
    setSubParticleIds(std::next(unordered.begin(), startSize), unordered.end());
  }

  ACTS_DEBUG("Found " << unordered.size() << " particles");

  // We need to keep the original indices, so that we can set the parent/child
  // relationship later
  std::vector<std::size_t> indices(unordered.size());
  std::iota(indices.begin(), indices.end(), 0);

  // Sort indices
  std::sort(
      indices.begin(), indices.end(), [&](std::size_t lhs, std::size_t rhs) {
        return detail::CompareParticleId{}(unordered[lhs], unordered[rhs]);
      });

  SimParticleContainer::sequence_type ordered;
  ordered.reserve(unordered.size());

  std::vector<std::size_t> unorderedToOrdered(indices.size());
  std::size_t unorderedIdx = 0;
  for (std::size_t idx : indices) {
    ordered.push_back(unordered[idx]);
    unorderedToOrdered[idx] = unorderedIdx;
    unorderedIdx++;
  }

  // ParentRelationship orderedParentRelationship;

  // Write ordered particles container to the EventStore
  SimParticleContainer particles;
  // This should be O(N) since the input is already ordered
  particles.insert(ordered.begin(), ordered.end());

  // Now that addresses should be stable, set parent/child relationships from
  // the indices we have above
  // Indices are in the unordered container, so need to convert to ordered
  for (const auto [unorderedChildIdx, unorderedParentIdx] :
       parentRelationship) {
    std::size_t orderedChildIdx = unorderedToOrdered[unorderedChildIdx];
    std::size_t orderedParentIdx = unorderedToOrdered[unorderedParentIdx];

    // get the pointers from the final ordered container
    auto& parent = *particles.find(ordered[orderedParentIdx].particleId());
    auto& child = *particles.find(ordered[orderedChildIdx].particleId());

    parent.addChild(child);
    child.setParent(&parent);
  }

  std::ofstream dot("particles.dot");
  graphvizSimParticleContainer(dot, particles);

  m_outputParticles(ctx, std::move(particles));

  return ProcessCode::SUCCESS;
}

void EDM4hepReader::processChildren(
    const edm4hep::MCParticle& inParticle, SimBarcode parentId,
    SimParticleContainer::sequence_type& particles,
    ParentRelationship& parentRelationship, std::size_t& nSecondaryVertices,
    std::size_t& maxGen) const {
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
    // std::cout << "DOT: "
    // << "P_" << inParticle.id().index << " -- "
    // << "P_" << daughter.id().index << "[label=\"" << particle.pdg()
    // << "\"]" << std::endl;

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
    parentRelationship[particles.size() - 1] = parentIndex;
    processChildren(daughter, pid, particles, parentRelationship,
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
