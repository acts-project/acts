// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <string>

#include "ACTFW/Framework/BareAlgorithm.hpp"

namespace FW {

/// Print hits within some geometric region-of-interest.
class PrintHits : public BareAlgorithm {
 public:
  struct Config {
    /// Input cluster collection.
    std::string inputClusters;
    /// Input hit-particles map.
    std::string inputHitParticlesMap;
    /// Input hit id collection
    std::string inputHitIds;
    // the following values are (invalid) defaults; should be set by user
    /// Hit id range which should be printed.
    size_t hitIdStart = 0u;
    size_t hitIdLength = 0u;
    /// Detector selection for which input will be printed.
    size_t volumeId = 0u;
    size_t layerId = 0u;
    size_t moduleId = 0u;
  };

  PrintHits(const Config& cfg, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const;

 private:
  Config m_cfg;
};

}  // namespace FW
