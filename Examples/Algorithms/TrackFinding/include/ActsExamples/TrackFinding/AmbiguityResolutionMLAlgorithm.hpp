// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Onnx/OnnxRuntimeBase.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

/// Evicts tracks that seem to be duplicated and fake.
///
/// The implementation works as follows:
///  1) Cluster together nearby tracks using shared hits
///  2) For each track use a neural network to compute a score
///  3) In each cluster keep the track with the highest score
class AmbiguityResolutionMLAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input trajectories collection.
    std::string inputTrajectories;
    /// path to the ONNX model for the duplicate neural network
    std::string inputDuplicateNN;
    /// Output trajectories collection.
    std::string outputTrajectories;
    /// Minumum number of measurement to form a track.
    int nMeasurementsMin = 7;
  };

  /// Construct the ambiguity resolution algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  AmbiguityResolutionMLAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Run the ambiguity resolution algorithm.
  ///
  /// @param cxt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  // ONNX environement
  Ort::Env m_env;
  // ONNX model for the duplicate neural network
  Acts::OnnxRuntimeBase m_duplicateClassifier;
};

}  // namespace ActsExamples
