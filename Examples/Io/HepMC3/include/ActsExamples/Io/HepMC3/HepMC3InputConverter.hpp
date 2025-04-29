// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace HepMC3 {
class GenEvent;
class GenVertex;
}  // namespace HepMC3

namespace ActsExamples {

class HepMC3InputConverter : public IAlgorithm {
 public:
  struct Config {
    std::string inputEvent;
    std::string outputParticles;
    std::string outputVertices;

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
  };

  HepMC3InputConverter(const Config& config, Acts::Logging::Level level);

  const Config& config() const { return m_cfg; }

  ProcessCode execute(const AlgorithmContext& ctx) const final;

 private:
  void convertHepMC3ToInternalEdm(const AlgorithmContext& ctx,
                                  const HepMC3::GenEvent& genEvent) const;

  void handleVertex(const HepMC3::GenVertex& genVertex, SimVertex& vertex,
                    std::vector<SimVertex>& vertices,
                    std::vector<SimParticle>& particles,
                    std::size_t& nSecondaryVertices, std::size_t& nParticles,
                    std::vector<bool>& seenVertices) const;

  Config m_cfg;

  ReadDataHandle<std::shared_ptr<HepMC3::GenEvent>> m_inputEvent{this,
                                                                 "InputEvent"};

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
  WriteDataHandle<SimVertexContainer> m_outputVertices{this, "OutputVertices"};
};

}  // namespace ActsExamples
