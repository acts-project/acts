// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/AmbiguityResolutionAlgorithm.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iterator>
#include <numeric>
#include <stdexcept>

#include <boost/container/flat_set.hpp>

ActsExamples::AmbiguityResolutionAlgorithm::AmbiguityResolutionAlgorithm(
    ActsExamples::AmbiguityResolutionAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("AmbiguityResolutionAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
}

namespace {

struct State {
  std::vector<std::pair<size_t, size_t>> trackTips;
  std::vector<float> trackChi2;
  std::vector<ActsExamples::TrackParameters> trackParameters;
  std::vector<std::vector<std::size_t>> measurementsPerTrack;

  std::vector<boost::container::flat_set<std::size_t>> tracksPerMeasurement;
  std::vector<std::size_t> sharedMeasurementsPerTrack;

  boost::container::flat_set<std::size_t> selectedTracks;
};

State computeInitialState(
    const ActsExamples::TrajectoriesContainer& trajectories,
    std::size_t nMeasurementsMin) {
  State state;

  {
    std::size_t iTrack = 0;
    for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
      const auto& traj = trajectories[iTraj];
      for (auto tip : traj.tips()) {
        if (!traj.hasTrackParameters(tip)) {
          continue;
        }

        auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
            traj.multiTrajectory(), tip);
        if (trajState.nMeasurements < nMeasurementsMin) {
          continue;
        }

        std::vector<std::size_t> measurements;
        traj.multiTrajectory().visitBackwards(tip, [&](const auto& hit) {
          if (hit.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
            std::size_t iMeasurement =
                hit.getUncalibratedSourceLink()
                    .template get<ActsExamples::IndexSourceLink>()
                    .index();
            measurements.push_back(iMeasurement);
          }
          return true;
        });

        state.trackTips.emplace_back(iTraj, tip);
        state.trackChi2.push_back(trajState.chi2Sum / trajState.NDF);
        state.trackParameters.push_back(traj.trackParameters(tip));
        state.measurementsPerTrack.push_back(std::move(measurements));

        state.selectedTracks.insert(iTrack);

        ++iTrack;
      }
    }
  }

  state.measurementsPerTrack = std::vector<std::vector<std::size_t>>(
      state.trackTips.size(), std::vector<std::size_t>());

  for (std::size_t iTrack = 0; iTrack < state.trackTips.size(); ++iTrack) {
    for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
      state.tracksPerMeasurement[iMeasurement].insert(iTrack);
    }
  }

  state.sharedMeasurementsPerTrack =
      std::vector<std::size_t>(state.trackTips.size(), 0);

  for (std::size_t iTrack = 0; iTrack < state.trackTips.size(); ++iTrack) {
    for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
      if (state.tracksPerMeasurement[iMeasurement].size() > 1) {
        ++state.sharedMeasurementsPerTrack[iTrack];
      }
    }
  }

  return state;
}

void removeTrack(State& state, std::size_t iTrack) {
  for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
    state.tracksPerMeasurement[iMeasurement].erase(iTrack);

    if (state.tracksPerMeasurement[iMeasurement].size() <= 1) {
      --state.sharedMeasurementsPerTrack[iTrack];
    }
  }

  state.selectedTracks.erase(iTrack);
}

}  // namespace

ActsExamples::ProcessCode ActsExamples::AmbiguityResolutionAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  auto state = computeInitialState(trajectories, m_cfg.nMeasurementsMin);

  auto sharedMeasurementsComperator = [&state](std::size_t a, std::size_t b) {
    return state.sharedMeasurementsPerTrack[a] -
           state.sharedMeasurementsPerTrack[b];
  };
  auto badTrackComperator = [&state](std::size_t a, std::size_t b) {
    auto loss = [&state](std::size_t i) {
      return 1.0f * state.sharedMeasurementsPerTrack[i] /
             state.measurementsPerTrack[i].size();
    };
    return loss(a) - loss(b);
  };

  while (true) {
    if (*std::max_element(
            std::begin(state.selectedTracks), std::end(state.selectedTracks),
            sharedMeasurementsComperator) < m_cfg.maximumSharedHits) {
      break;
    }

    const auto badTrack =
        *std::max_element(std::begin(state.selectedTracks),
                          std::end(state.selectedTracks), badTrackComperator);
    removeTrack(state, badTrack);
  }

  ACTS_INFO("Resolved to " << state.selectedTracks.size() << " tracks from "
                           << state.trackTips.size());

  TrajectoriesContainer outputTrajectories;
  outputTrajectories.reserve(trajectories.size());
  for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
    const auto& traj = trajectories[iTraj];

    std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;
    Trajectories::IndexedParameters parameters;

    for (auto iTrack : state.selectedTracks) {
      if (state.trackTips[iTrack].first != iTraj) {
        continue;
      }
      const auto tip = state.trackTips[iTrack].second;
      tips.push_back(tip);
      parameters.emplace(tip, state.trackParameters[iTrack]);
    }
    if (!tips.empty()) {
      outputTrajectories.emplace_back(traj.multiTrajectory(), tips, parameters);
    }
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(outputTrajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}
