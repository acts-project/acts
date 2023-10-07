// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/TrajectoriesToPrototracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

namespace ActsExamples {

TrajectoriesToPrototracks::TrajectoriesToPrototracks(Config cfg,
                                                     Acts::Logging::Level lvl)
    : IAlgorithm("TrajectoriesToPrototracks", lvl), m_cfg(std::move(cfg)) {
  m_inputTrajectories.initialize(m_cfg.inputTrajectories);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode TrajectoriesToPrototracks::execute(
    const AlgorithmContext& ctx) const {
  const auto trajectories = m_inputTrajectories(ctx);

  ProtoTrackContainer tracks;

  for (const auto& trajectory : trajectories) {
    for (const auto tip : trajectory.tips()) {
      ProtoTrack track;

      trajectory.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
        if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        auto source_link =
            state.getUncalibratedSourceLink().template get<IndexSourceLink>();
        track.push_back(source_link.index());

        return true;
      });

      tracks.push_back(track);
    }
  }

  m_outputProtoTracks(ctx, std::move(tracks));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
