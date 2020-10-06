// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>

namespace ActsExamples {

/// Print all particles.
class ParticlesPrinter : public BareAlgorithm {
 public:
  struct Config {
    /// Input particles collection.
    std::string inputParticles;
  };

  ParticlesPrinter(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
