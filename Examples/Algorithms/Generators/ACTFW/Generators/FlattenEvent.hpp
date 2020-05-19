// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// Convert nested vector of vertices into a vector of particles.
///
/// Ignores all incoming particles in the vertices.
class FlattenEvent final : public BareAlgorithm {
 public:
  struct Config {
    /// The input event collection.
    std::string inputEvent;
    /// The output particles collection.
    std::string outputParticles;
  };

  FlattenEvent(const Config& cfg, Acts::Logging::Level lvl);

  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
