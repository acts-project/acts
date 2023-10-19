// This file is part of the Acts project.
//
// Copyright (C) 2017-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Required.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>

namespace ActsExamples {

/// Generate particles from uniform parameter distributions.
///
/// Generates a single vertex with the given number of tracks. The track
/// direction is drawn from a uniform distribution on the unit sphere (within
/// the given limits). Its absolute momentum is drawn from a uniform
/// distribution. Position and time are always set to zero.
class ParametricParticleGenerator : public EventGenerator::ParticlesGenerator {
 public:
  struct Config {
    /// Low, high (exclusive) for the transverse direction angle.
    double phiMin = -M_PI;
    double phiMax = M_PI;

    /// Low, high (exclusive) for the pseudo rapidity.
    ///
    /// @note This is the standard generation, for detector performance
    /// classification, where a flat distribution in eta can be useful,
    /// this can be set by the etaUniform flag;
    ///
    Required<double> etaMin;
    Required<double> etaMax;
    Required<bool> etaUniform;

    /// Low, high (exclusive) for absolute/transverse momentum.
    Required<double> pMin;
    Required<double> pMax;
    /// Indicate if the momentum referse to transverse momentum.
    Required<bool> pTransverse;

    /// (Absolute) PDG particle number to identify the particle type.
    Required<Acts::PdgParticle> pdg;
    /// Randomize the charge and flip the PDG particle number sign accordingly.
    Required<bool> randomizeCharge;

    /// Number of particles.
    size_t numParticles = 1;

    /// Overrides particle charge.
    std::optional<double> charge;
    /// Overrides particle mass.
    std::optional<double> mass;
  };

  ParametricParticleGenerator(const Config& cfg);

  /// Generate a single primary vertex with the given number of particles.
  SimParticleContainer operator()(RandomEngine& rng) override;

 private:
  Config m_cfg;

  double m_cosThetaMin{};
  double m_cosThetaMax{};

  double m_charge{};
  double m_mass{};
};

}  // namespace ActsExamples
