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
#include "ActsExamples/TruthTracking/ParticleSmearingOptions.hpp"

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
    void ReadOptions(const Options::Variables& vars) {
      using namespace Acts::UnitConstants;
      using Options::Reals;

      auto sigmaD0Opts = vars["smear-sigma-D0"].template as<Reals<3>>();
      auto sigmaZ0Opts = vars["smear-sigma-Z0"].template as<Reals<3>>();
      auto sigmaMomOpts = vars["smear-sigma-momentum"].template as<Reals<3>>();

      sigmaD0 = sigmaD0Opts[0] * Acts::UnitConstants::um;
      sigmaD0PtA = sigmaD0Opts[1] * Acts::UnitConstants::um;
      sigmaD0PtB = sigmaD0Opts[2] / Acts::UnitConstants::GeV;
      sigmaZ0 = sigmaZ0Opts[0] * Acts::UnitConstants::um;
      sigmaZ0PtA = sigmaZ0Opts[1] * Acts::UnitConstants::um;
      sigmaZ0PtB = sigmaZ0Opts[2] / Acts::UnitConstants::GeV;
      sigmaT0 = vars["smear-sigma-T0"].as<double>() * Acts::UnitConstants::ns;
      sigmaPhi = sigmaMomOpts[0] * Acts::UnitConstants::degree;
      sigmaTheta = sigmaMomOpts[1] * Acts::UnitConstants::degree;
      sigmaPRel = sigmaMomOpts[2];
      initialVarInflation =
          vars["fit-initial-variance-inflation"].template as<Reals<6>>();
    }

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
    /// Inflate the initial covariance matrix
    std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
    /// Random numbers service.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
  };

  ParticleSmearing(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
