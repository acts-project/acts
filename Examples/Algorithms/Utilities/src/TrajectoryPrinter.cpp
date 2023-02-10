// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Utilities/TrajectoryPrinter.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

namespace ActsExamples {

ProcessCode TrajectoryPrinter::execute(const AlgorithmContext& ctx) const {
  const auto trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  std::vector<unsigned int> nTracksPerSeedCount;
  nTracksPerSeedCount.push_back(0);

  for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
    ACTS_INFO("- trajectory #" << iTraj);
    const auto& traj = trajectories[iTraj];
    const auto& mj = traj.multiTrajectory();

    if (traj.tips().size() + 1 >= nTracksPerSeedCount.size()) {
      nTracksPerSeedCount.resize(traj.tips().size() + 1);
    }
    nTracksPerSeedCount.at(traj.tips().size()) += 1;

    ACTS_INFO(" -> " << traj.tips().size() << " tips");
    for (auto trackTip : traj.tips()) {
      const auto& fittedParameters = traj.trackParameters(trackTip);
      ACTS_INFO(" --> tip: " << trackTip << ""
                             << " with "
                             << fittedParameters.parameters().transpose());
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

      ACTS_INFO(" --> tip state:");
      ACTS_INFO("     - nStates: " << trajState.nStates);
      ACTS_INFO("     - nMeasurements: " << trajState.nMeasurements);
      ACTS_INFO("     - nOutliers: " << trajState.nOutliers);
      ACTS_INFO("     - nHoles: " << trajState.nHoles);
      ACTS_INFO("     - chi2Sum: " << trajState.chi2Sum);
      ACTS_INFO("     - NDF: " << trajState.NDF);
      ACTS_INFO("     - nSharedHits: " << trajState.nSharedHits);

      // for (const auto ts : mj.trackStateRange(trackTip)) {
      // }
    }
  }

  ACTS_INFO("" << trajectories.size() << " trajectories total");
  for (size_t idx = 0; idx < nTracksPerSeedCount.size(); idx++) {
    ACTS_INFO(" - " << idx << " tracks -> " << nTracksPerSeedCount.at(idx)
                    << " times");
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
