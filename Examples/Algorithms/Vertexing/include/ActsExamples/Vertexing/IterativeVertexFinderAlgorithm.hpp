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
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

class IterativeVertexFinderAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Optional. Input track parameters collection
    std::string inputTrackParameters;
    /// Optional. Input trajectories container.
    std::string inputTrajectories;
    /// Output proto vertex collection
    std::string outputProtoVertices;
    /// Output vertex collection
    std::string outputVertices = "vertices";
    /// Output reconstruction time in ms
    std::string outputTime = "time";
    /// The magnetic field
    std::shared_ptr<Acts::MagneticFieldProvider> bField;
  };

  IterativeVertexFinderAlgorithm(const Config& config,
                                 Acts::Logging::Level level);

  ~IterativeVertexFinderAlgorithm();

  /// Find vertices using iterative vertex finder algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  ProcessCode initialize() final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
