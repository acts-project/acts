// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/TrajectoriesToProtoTracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include <utility>
#include <vector>

namespace ActsExamples {

TrajectoriesToProtoTracks::TrajectoriesToProtoTracks(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("TrajectoriesToProtoTracks", std::move(logger)),
      m_cfg(std::move(cfg)) {
  m_inputTrajectories.initialize(m_cfg.inputTrajectories);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode TrajectoriesToProtoTracks::execute(
    const AlgorithmContext& ctx) const {
  const auto trajectories = m_inputTrajectories(ctx);

  ProtoTrackContainer tracks;

  for (const auto& trajectory : trajectories) {
    for (const auto tip : trajectory.tips()) {
      ProtoTrack track;

      trajectory.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
        if (!state.typeFlags().isMeasurement()) {
          return true;
        }

        const auto sourceLink =
            state.getUncalibratedSourceLink().template get<IndexSourceLink>();
        track.push_back(sourceLink.index());

        return true;
      });

      tracks.push_back(track);
    }
  }

  m_outputProtoTracks(ctx, std::move(tracks));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
