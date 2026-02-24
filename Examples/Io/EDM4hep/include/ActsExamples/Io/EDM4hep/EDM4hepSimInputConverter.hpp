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
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"
#include "ActsPlugins/EDM4hep/EDM4hepUtil.hpp"

#include <memory>
#include <string>

namespace edm4hep {
class MCParticle;
class SimTrackerHit;
}  // namespace edm4hep

namespace podio {
class Frame;
}

namespace ActsExamples {

namespace detail {
struct ParticleInfo;
}

class DD4hepDetector;

using EDM4hepSimHitAssociation = std::vector<edm4hep::SimTrackerHit>;

/// Read particles from EDM4hep.
///
/// Inpersistent information:
/// - particle ID
/// - process
class EDM4hepSimInputConverter final : public PodioInputConverter {
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
    /// Output a mapping from internal sim hit index to edm4hep input hits
    /// @note Optional, will not be computed if it's not stored
    std::optional<std::string> outputSimHitAssociation = std::nullopt;
    /// Output simulated vertices collection.
    std::string outputSimVertices;

    /// DD4hep detector for cellID resolution.
    std::shared_ptr<DD4hepDetector> dd4hepDetector;

    /// Tracking geometry for cellID resolution.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// Whether to sort sim hits in time to produce index sequence
    bool sortSimHitsInTime = false;

    std::optional<double> particleRMin = std::nullopt;
    std::optional<double> particleRMax = std::nullopt;

    std::optional<double> particleZMin = std::nullopt;
    std::optional<double> particleZMax = std::nullopt;

    std::optional<double> particlePtMin = std::nullopt;
    std::optional<double> particlePtMax = std::nullopt;
  };

  using ParentRelationship = std::unordered_map<std::size_t, std::size_t>;

  /// Construct the particle reader.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  explicit EDM4hepSimInputConverter(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode convert(const AlgorithmContext& ctx,
                      const podio::Frame& frame) const override;

  void processChildren(
      const edm4hep::MCParticle& particle, SimBarcode parentId,
      std::vector<SimBarcode>& particles,
      ParentRelationship& parentRelationship,
      std::unordered_map<int, detail::ParticleInfo>& particleMap,
      std::size_t& nSecondaryVertices, std::size_t& maxGen,
      const std::function<std::uint8_t(const edm4hep::MCParticle&)>& getNumHits)
      const;

  static void setSubParticleIds(std::span<SimBarcode> particles);

  bool acceptParticle(const edm4hep::MCParticle& particle) const;

  bool particleOrDescendantsHaveHits(
      const edm4hep::MCParticle& particle,
      const std::function<std::uint8_t(const edm4hep::MCParticle&)>& getNumHits)
      const;

  Config m_cfg;

  std::unordered_map<unsigned int, const Acts::Surface*> m_surfaceMap;

  WriteDataHandle<SimParticleContainer> m_outputParticlesGenerator{
      this, "OutputParticlesGenerator"};
  WriteDataHandle<SimParticleContainer> m_outputParticlesSimulation{
      this, "OutputParticlesSimulation"};

  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHits"};
  WriteDataHandle<ActsPlugins::EDM4hepUtil::SimHitAssociation>
      m_outputSimHitAssociation{this, "OutputSimHitAssociation"};

  WriteDataHandle<SimVertexContainer> m_outputSimVertices{this,
                                                          "OutputSimVertices"};
};

}  // namespace ActsExamples
