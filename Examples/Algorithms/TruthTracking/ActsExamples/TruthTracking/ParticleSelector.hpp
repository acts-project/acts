// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @date 2018-03-14
/// @author Moritz Kiehn <msmk@cern.ch>

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <limits>

namespace ActsExamples {

/// Select particles by applying some selection cuts.
class ParticleSelector final : public BareAlgorithm {
 public:
  struct Config {
    /// The input particles collection.
    std::string inputParticles;
    /// The output particles collection.
    std::string outputParticles;
    // Minimum/maximum distance from the origin in the tranverse plane.
    double rhoMin = 0;
    double rhoMax = std::numeric_limits<double>::infinity();
    // Minimum/maximum absolute distance from the origin along z.
    double absZMin = 0;
    double absZMax = std::numeric_limits<double>::infinity();
    // Minimum/maximum particle time.
    double timeMin = -std::numeric_limits<double>::infinity();
    double timeMax = std::numeric_limits<double>::infinity();
    // Direction cuts.
    double phiMin = -std::numeric_limits<double>::infinity();
    double phiMax = std::numeric_limits<double>::infinity();
    double etaMin = -std::numeric_limits<double>::infinity();
    double etaMax = std::numeric_limits<double>::infinity();
    double absEtaMin = 0;
    double absEtaMax = std::numeric_limits<double>::infinity();
    // Momentum cuts.
    double ptMin = 0;
    double ptMax = std::numeric_limits<double>::infinity();
    /// Remove charged particles.
    bool removeCharged = false;
    /// Remove neutral particles.
    bool removeNeutral = false;
  };

  ParticleSelector(const Config& config, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
