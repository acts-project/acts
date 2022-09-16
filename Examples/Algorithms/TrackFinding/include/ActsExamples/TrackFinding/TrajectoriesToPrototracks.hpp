// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/BareAlgorithm.hpp"

namespace ActsExamples {

class TrajectoriesToPrototracks final : public BareAlgorithm {
 public:
  struct Config {
    std::string inputTrajectories = "trajectories";
    std::string outputPrototracks = "tracks-from-trajectories";
  };

  /// Construct the algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  TrajectoriesToPrototracks(Config cfg, Acts::Logging::Level lvl)
      : BareAlgorithm("TrajectoriesToPrototracks", lvl), m_cfg(cfg) {}

  /// Run the algorithm.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
