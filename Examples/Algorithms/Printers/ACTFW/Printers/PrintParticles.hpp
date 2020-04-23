// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// Print all particles.
class PrintParticles : public BareAlgorithm {
 public:
  struct Config {
    /// Input particles collection.
    std::string inputParticles;
  };

  PrintParticles(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const;

 private:
  Config m_cfg;
};

}  // namespace FW
