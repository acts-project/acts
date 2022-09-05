// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/TrajectoriesToPrototracks.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

namespace ActsExamples {

ProcessCode TrajectoriesToPrototracks::execute(
    const AlgorithmContext& ctx) const {
  const auto trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  ProtoTrackContainer tracks;

  for (const auto& trajectory : trajectories) {
    for (const auto tip : trajectory.tips()) {
      ProtoTrack track;

      trajectory.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
        if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        const auto& source_link =
            static_cast<const IndexSourceLink&>(state.uncalibrated());
        track.push_back(source_link.index());

        return true;
      });

      tracks.push_back(track);
    }
  }

  ctx.eventStore.add(m_cfg.outputPrototracks, std::move(tracks));

  return ProcessCode::SUCCESS;
}
}  // namespace ActsExamples
