// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <cstddef>
#include <string>

namespace ActsExamples {

/// Print hits within some geometric region-of-interest.
class HitsPrinter : public BareAlgorithm {
 public:
  struct Config {
    /// Input cluster collection.
    std::string inputClusters;
    /// Input hit-particles map.
    std::string inputMeasurementParticlesMap;
    /// Input hit id collection
    std::string inputHitIds;
    // Print hits selected by their indices (zero length to disable).
    size_t selectIndexStart = 0u;
    size_t selectIndexLength = 0u;
    // Print hits within a certain geometry range (zero to disable).
    size_t selectVolume = 0u;
    size_t selectLayer = 0u;
    size_t selectModule = 0u;
  };

  HitsPrinter(const Config& cfg, Acts::Logging::Level level);

  ProcessCode execute(const AlgorithmContext& ctx) const;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
