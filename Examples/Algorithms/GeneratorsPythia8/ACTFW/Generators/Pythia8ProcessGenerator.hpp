// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <mutex>
#include <vector>

#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Pythia8 {
class Pythia;
}
namespace FW {

class Pythia8Generator {
 public:
  struct Config {
    /// PDG particle number of the first incoming beam.
    Acts::PdgParticle pdgBeam0 = Acts::PdgParticle::eProton;
    /// PDG particle number of the second incoming beam.
    Acts::PdgParticle pdgBeam1 = Acts::PdgParticle::eProton;
    /// Center-of-mass energy.
    double cmsEnergy = 14 * Acts::UnitConstants::TeV;
    /// Additional Pythia8 settings.
    std::vector<std::string> settings = {{"HardQCD:all = on"}};
  };

  static std::function<std::vector<SimVertex>(RandomEngine&)> makeFunction(
      const Config& cfg, Acts::Logging::Level lvl);

  // try to prevent pythia breakage by forbidding copying
  Pythia8Generator() = delete;
  Pythia8Generator(const Pythia8Generator&) = delete;
  Pythia8Generator(Pythia8Generator&&) = delete;
  Pythia8Generator& operator=(const Pythia8Generator&) = delete;
  Pythia8Generator& operator=(Pythia8Generator&& other) = delete;

  Pythia8Generator(const Config& cfg, Acts::Logging::Level lvl);
  ~Pythia8Generator();

  std::vector<SimVertex> operator()(RandomEngine& rng);

 private:
  /// Private access to the logging instance
  const Acts::Logger& logger() const { return (*m_logger); }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<::Pythia8::Pythia> m_pythia8;
  std::mutex m_pythia8Mutex;
};

}  // namespace FW
