// This file is part of the Acts project.
//
// Copyright (C) 2016-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>

namespace ActsExamples {

class IterativeVertexFinderAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input track parameters collection
    std::string inputTrackParameters;
    /// Output proto vertex collection
    std::string outputProtoVertices;
    /// Output vertex collection
    std::string outputVertices = "vertices";
    /// Output reconstruction time in ms
    std::string outputTime = "time";
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
  };

  IterativeVertexFinderAlgorithm(const Config& cfg, Acts::Logging::Level level);

  /// Find vertices using iterative vertex finder algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
