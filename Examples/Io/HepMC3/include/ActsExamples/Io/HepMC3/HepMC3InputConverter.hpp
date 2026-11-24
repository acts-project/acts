// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticleDecay.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <optional>
#include <string>
#include <vector>

namespace HepMC3 {
class GenEvent;
class GenVertex;
class GenParticle;
}  // namespace HepMC3

namespace ActsExamples {

class HepMC3InputConverter : public IAlgorithm {
 public:
  struct Config {
    std::string inputEvent;
    std::string outputParticles;
    std::string outputVertices;

    /// Optional decay graph; outputParticles then holds only transport roots.
    std::string outputDecayGraph;
    /// Signed PDG IDs whose physical HepMC decays should be prescribed.
    std::vector<int> decayPdgCodes;

    /// If true, print the listing of the generated event. This can be very
    /// verbose
    bool printListing = false;

    /// Merge primary vertices
    bool mergePrimaries = true;

    /// The spatial vertex threshold below which to consider primary vertices
    /// candidates identical.
    double primaryVertexSpatialThreshold = 1 * Acts::UnitConstants::nm;

    /// The spatial vertex threshold below which to consider secondary vertices
    /// candidates identical.
    double vertexSpatialThreshold = 1 * Acts::UnitConstants::um;

    /// If true, merge secondary vertices that are close to their parent vertex
    bool mergeSecondaries = true;

    /// If true, check the consistency of the generated event.
    bool checkConsistency = false;

    /// Search up to quark in heavy-flavour tagging
    bool searchUpToHeavyFlavourQuark = false;
  };

  explicit HepMC3InputConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  const Config& config() const { return m_cfg; }

  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  void convertHepMC3ToInternalEdm(const AlgorithmContext& ctx,
                                  const HepMC3::GenEvent& genEvent) const;

  HeavyFlavourOrigin deriveHeavyFlavourOrigin(
      const std::shared_ptr<const HepMC3::GenParticle>& particleToCheck) const;

  void handleVertex(const HepMC3::GenVertex& genVertex, SimVertex& vertex,
                    std::vector<SimVertex>& vertices,
                    std::vector<SimParticle>& particles,
                    std::size_t& nSecondaryVertices, std::size_t& nParticles,
                    std::vector<bool>& seenVertices,
                    SimParticleDecayGraph* graph) const;

  std::optional<SimParticle> makeParticle(
      const std::shared_ptr<const HepMC3::GenParticle>& source, SimBarcode id,
      const Acts::Vector4& position, bool prescribed) const;

  void addDecay(const HepMC3::GenParticle& source, const SimParticle& parent,
                SimParticleDecayGraph& graph, std::vector<bool>& seenVertices,
                SimBarcode::SubParticleId& nextSubparticle) const;

  Config m_cfg;

  ReadDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_inputEvent{this,
                                                                 "InputEvent"};

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
  WriteDataHandle<SimParticleDecayGraph> m_outputDecayGraph{this,
                                                            "OutputDecayGraph"};
  WriteDataHandle<SimVertexContainer> m_outputVertices{this, "OutputVertices"};
};

}  // namespace ActsExamples
