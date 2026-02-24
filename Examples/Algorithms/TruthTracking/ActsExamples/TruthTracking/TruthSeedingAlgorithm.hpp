// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ActsExamples {

/// Construct track seeds from particles.
class TruthSeedingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// The input truth particles that should be used for truth seeding.
    std::string inputParticles;
    /// The input particle-measurements map collection.
    std::string inputParticleMeasurementsMap;
    /// The input sim hits collection that is used to order space points in the
    /// seeds.
    std::string inputSimHits;
    /// The input measurement-sim hits map collection.
    std::string inputMeasurementSimHitsMap;
    /// Input space point collections.
    ///
    /// We allow multiple space point collections to allow different parts of
    /// the detector to use different algorithms for space point construction,
    /// e.g. single-hit space points for pixel-like detectors or double-hit
    /// space points for strip-like detectors.
    std::vector<std::string> inputSpacePoints;
    /// Output successfully seeded truth particles.
    std::string outputParticles;
    /// Output proto track collection.
    std::string outputProtoTracks;
    /// Output seed collection.
    std::string outputSeeds;
    /// Optional. Output particle hypotheses collection.
    std::string outputParticleHypotheses;

    /// Optional particle hypothesis override.
    std::optional<Acts::ParticleHypothesis> particleHypothesis = std::nullopt;

    /// Minimum deltaR between space points in a seed
    float deltaRMin = 10 * Acts::UnitConstants::mm;
    /// Maximum deltaR between space points in a seed
    float deltaRMax = 200 * Acts::UnitConstants::mm;
    /// Minimum absDeltaZMin between space points in a seed
    float absDeltaZMin = 0 * Acts::UnitConstants::mm;
    /// Maximum absDeltaZMax between space points in a seed
    float absDeltaZMax = 500 * Acts::UnitConstants::mm;
  };

  /// Construct the truth seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  explicit TruthSeedingAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Run the truth seeding algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<InverseMultimap<SimBarcode>> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputHits"};
  ReadDataHandle<InverseMultimap<Index>> m_inputMeasurementSimHitsMap{
      this, "MeasurementSimHitsMap"};
  std::vector<std::unique_ptr<ReadDataHandle<SimSpacePointContainer>>>
      m_inputSpacePoints{};

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
  WriteDataHandle<ProtoTrackContainer> m_outputProtoTracks{this,
                                                           "OutputProtoTracks"};
  WriteDataHandle<SimSeedContainer> m_outputSeeds{this, "OutputSeeds"};
  WriteDataHandle<std::vector<Acts::ParticleHypothesis>>
      m_outputParticleHypotheses{this, "OutputParticleHypotheses"};
};

}  // namespace ActsExamples
