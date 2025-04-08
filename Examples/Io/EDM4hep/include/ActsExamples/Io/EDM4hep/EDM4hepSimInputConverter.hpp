// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepInputConverter.hpp"

#include <memory>
#include <string>

namespace edm4hep {
class MCParticle;
}

namespace podio {
class Frame;
}

namespace ActsExamples {

class DD4hepDetector;

/// Read particles from EDM4hep.
///
/// Inpersistent information:
/// - particle ID
/// - process
class EDM4hepSimInputConverter final : public EDM4hepInputConverter {
 public:
  struct Config {
    /// Where to read input file from.
    std::string inputFrame = "events";
    /// Name of the particle collection in EDM4hep.
    std::string inputParticles = "MCParticles";
    /// Names of the sim hit collections
    std::vector<std::string> inputSimHits{};
    /// Particles from the generator
    std::string outputParticlesGenerator;
    /// Particles from the simulation
    std::string outputParticlesSimulation;
    /// Output simulated (truth) hits collection.
    std::string outputSimHits;
    /// Output simulated vertices collection.
    std::string outputSimVertices;

    /// Directory into which to write graphviz files for particles
    /// Empty string means no output
    std::string graphvizOutput = "";

    /// DD4hep detector for cellID resolution.
    std::shared_ptr<DD4hepDetector> dd4hepDetector;

    /// Tracking geometry for cellID resolution.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Whether to sort sim hits in time to produce index sequence
    bool sortSimHitsInTime = false;
  };

  using ParentRelationship = std::unordered_map<std::size_t, std::size_t>;

  /// Construct the particle reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  EDM4hepSimInputConverter(const Config& config, Acts::Logging::Level level);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const override;

  void processChildren(const edm4hep::MCParticle& particle, SimBarcode parentId,
                       std::vector<SimParticle>& particles,
                       ParentRelationship& parentRelationship,
                       std::unordered_map<int, std::size_t>& particleMap,
                       std::size_t& nSecondaryVertices,
                       std::size_t& maxGen) const;

  static void setSubParticleIds(std::span<SimParticle> particles);

  Config m_cfg;

  std::unordered_map<unsigned int, const Acts::Surface*> m_surfaceMap;

  WriteDataHandle<SimParticleContainer> m_outputParticlesGenerator{
      this, "OutputParticlesGenerator"};
  WriteDataHandle<SimParticleContainer> m_outputParticlesSimulation{
      this, "OutputParticlesSimulation"};

  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHits"};

  WriteDataHandle<SimVertexContainer> m_outputSimVertices{this,
                                                          "OutputSimVertices"};

  void graphviz(std::ostream& os, const std::vector<SimParticle>& particles,
                const ParentRelationship& parents) const;
};

}  // namespace ActsExamples
