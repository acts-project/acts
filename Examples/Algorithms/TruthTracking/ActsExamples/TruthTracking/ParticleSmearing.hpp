// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <string>

namespace ActsExamples {
class RandomNumbers;
struct AlgorithmContext;

/// Create track states by smearing truth particle information.
///
/// Particles are smeared in the perigee frame anchored at their true vertex
/// position. The `d0` and `z0` parameters are always defined within that
/// perigee frame and not globally. The generated bound parameters are stored in
/// the same order as the input particles.
class ParticleSmearing final : public IAlgorithm {
 public:
  struct Config {
    /// Input truth particles collection.
    std::string inputParticles;
    /// Output smeared tracks parameters collection.
    std::string outputTrackParameters;
    /// Constant term of the d0 resolution.
    double sigmaD0 = 20 * Acts::UnitConstants::um;
    /// Pt-dependent d0 resolution of the form sigma_d0 = A*exp(-1.*abs(B)*pt).
    double sigmaD0PtA = 30 * Acts::UnitConstants::um;
    double sigmaD0PtB = 0.3 / Acts::UnitConstants::GeV;
    /// Constant term of the z0 resolution.
    double sigmaZ0 = 20 * Acts::UnitConstants::um;
    /// Pt-dependent z0 resolution of the form sigma_z0 = A*exp(-1.*abs(B)*pt).
    double sigmaZ0PtA = 30 * Acts::UnitConstants::um;
    double sigmaZ0PtB = 0.3 / Acts::UnitConstants::GeV;
    /// Time resolution.
    double sigmaT0 = 1 * Acts::UnitConstants::ns;
    /// Phi angular resolution.
    double sigmaPhi = 1 * Acts::UnitConstants::degree;
    /// Theta angular resolution.
    double sigmaTheta = 1 * Acts::UnitConstants::degree;
    /// Relative momentum resolution.
    double sigmaPRel = 0.05;
    /// Optional. Initial covariance matrix diagonal. Overwrites the default if
    /// set.
    std::optional<std::array<double, 6>> initialSigmas = std::array<double, 6>{
        1 * Acts::UnitConstants::mm,     1 * Acts::UnitConstants::mm,
        1 * Acts::UnitConstants::degree, 1 * Acts::UnitConstants::degree,
        0.1 / Acts::UnitConstants::GeV,  1 * Acts::UnitConstants::ns};
    /// Inflate the initial covariance matrix
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
    /// Random numbers service.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
    /// Optional particle hypothesis override.
    std::optional<Acts::ParticleHypothesis> particleHypothesis = std::nullopt;
  };

  ParticleSmearing(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  WriteDataHandle<TrackParametersContainer> m_outputTrackParameters{
      this, "OutputTrackParameters"};
};

}  // namespace ActsExamples
