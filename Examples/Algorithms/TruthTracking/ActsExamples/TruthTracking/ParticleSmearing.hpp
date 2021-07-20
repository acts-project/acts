// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <array>
#include <limits>
#include <string>

namespace ActsExamples {

/// Create track states by smearing truth particle information.
///
/// Particles are smeared in the perigee frame anchored at their true vertex
/// position. The `d0` and `z0` parameters are always defined within that
/// perigee frame and not globally. The generated bound parameters are stored in
/// the same order as the input particles.
class ParticleSmearing final : public BareAlgorithm {
 public:
  struct Config {
    /// Input truth particles collection.
    std::string inputParticles;
    /// Output smeared tracks parameters collection.
    std::string outputTrackParameters;
    /// Constant term of the d0 resolution.
    double sigmaD0 = 30 * Acts::UnitConstants::um;
    /// Pt-dependent d0 resolution of the form sigma_d0 = A*exp(-1.*abs(B)*pt).
    double sigmaD0PtA = 0 * Acts::UnitConstants::um;
    double sigmaD0PtB = 1 / Acts::UnitConstants::GeV;
    /// Constant term of the z0 resolution.
    double sigmaZ0 = 30 * Acts::UnitConstants::um;
    /// Pt-dependent z0 resolution of the form sigma_z0 = A*exp(-1.*abs(B)*pt).
    double sigmaZ0PtA = 0 * Acts::UnitConstants::um;
    double sigmaZ0PtB = 1 / Acts::UnitConstants::GeV;
    /// Time resolution.
    double sigmaT0 = 5 * Acts::UnitConstants::ns;
    /// Phi angular resolution.
    double sigmaPhi = 1 * Acts::UnitConstants::degree;
    /// Theta angular resolution.
    double sigmaTheta = 1 * Acts::UnitConstants::degree;
    /// Relative momentum resolution.
    double sigmaPRel = 0.001;
    /// Inflate the initial covariance matrix
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
    /// Random numbers service.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
  };

  ParticleSmearing(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
