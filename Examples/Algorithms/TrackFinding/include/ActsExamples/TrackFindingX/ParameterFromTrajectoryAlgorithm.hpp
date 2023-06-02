// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

class ParameterFromTrajectoryAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input spacepoints collection.
    std::string inputTrajectories;

    /// Output protoTracks collection.
    std::string outputParamters;
  };

  /// Constructor of the track finding algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  ParameterFromTrajectoryAlgorithm(Config cfg, Acts::Logging::Level lvl)
      : BareAlgorithm("ParsFromTraj", lvl), m_cfg(cfg) {}

  virtual ~ParameterFromTrajectoryAlgorithm() {}

  /// Filter the measurements
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final {
    const auto& trajs =
        ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

    TrackParametersContainer trackParameters;

    for (const auto& traj : trajs) {
      const auto i = traj.tips().front();
      const auto state = traj.multiTrajectory().getTrackState(i);

      trackParameters.emplace_back(state.referenceSurface().getSharedPtr(),
                                   state.smoothed(),
                                   state.smoothedCovariance());
    }

    ctx.eventStore.add<TrackParametersContainer>(m_cfg.outputParamters,
                                                 std::move(trackParameters));

    return ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  // configuration
  Config m_cfg;
};

}  // namespace ActsExamples
