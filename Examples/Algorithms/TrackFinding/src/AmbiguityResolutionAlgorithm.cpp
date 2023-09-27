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

#include <boost/container/flat_map.hpp>
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

  m_inputTrajectories.initialize(m_cfg.inputTrajectories);
  m_outputTrajectories.initialize(m_cfg.outputTrajectories);
}

namespace {

struct State {
  std::size_t numberOfTracks{};

  std::vector<std::pair<std::size_t, std::size_t>> trackTips;
  std::vector<float> trackChi2;
  std::vector<ActsExamples::TrackParameters> trackParameters;
  std::vector<std::vector<std::size_t>> measurementsPerTrack;

  boost::container::flat_map<std::size_t,
                             boost::container::flat_set<std::size_t>>
      tracksPerMeasurement;
  std::vector<std::size_t> sharedMeasurementsPerTrack;

  boost::container::flat_set<std::size_t> selectedTracks;
};

State computeInitialState(
    const ActsExamples::TrajectoriesContainer& trajectories,
    std::size_t nMeasurementsMin) {
  State state;

  for (std::size_t iTrack = 0, iTraj = 0; iTraj < trajectories.size();
       ++iTraj) {
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

      ++state.numberOfTracks;

      state.trackTips.emplace_back(iTraj, tip);
      state.trackChi2.push_back(trajState.chi2Sum / trajState.NDF);
      state.trackParameters.push_back(traj.trackParameters(tip));
      state.measurementsPerTrack.push_back(std::move(measurements));

      state.selectedTracks.insert(iTrack);

      ++iTrack;
    }
  }

  for (std::size_t iTrack = 0; iTrack < state.numberOfTracks; ++iTrack) {
    for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
      state.tracksPerMeasurement[iMeasurement].insert(iTrack);
    }
  }

  state.sharedMeasurementsPerTrack =
      std::vector<std::size_t>(state.trackTips.size(), 0);

  for (std::size_t iTrack = 0; iTrack < state.numberOfTracks; ++iTrack) {
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

    if (state.tracksPerMeasurement[iMeasurement].size() == 1) {
      auto jTrack = *std::begin(state.tracksPerMeasurement[iMeasurement]);
      --state.sharedMeasurementsPerTrack[jTrack];
    }
  }

  state.selectedTracks.erase(iTrack);
}

}  // namespace

ActsExamples::ProcessCode ActsExamples::AmbiguityResolutionAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  const auto& trajectories = m_inputTrajectories(ctx);

  auto state = computeInitialState(trajectories, m_cfg.nMeasurementsMin);

  auto sharedMeasurementsComperator = [&state](std::size_t a, std::size_t b) {
    return state.sharedMeasurementsPerTrack[a] <
           state.sharedMeasurementsPerTrack[b];
  };
  auto badTrackComperator = [&state](std::size_t a, std::size_t b) {
    auto relativeSharedMeasurements = [&state](std::size_t i) {
      return 1.0 * state.sharedMeasurementsPerTrack[i] /
             state.measurementsPerTrack[i].size();
    };

    if (relativeSharedMeasurements(a) != relativeSharedMeasurements(b)) {
      return relativeSharedMeasurements(a) < relativeSharedMeasurements(b);
    }
    return state.trackChi2[a] < state.trackChi2[b];
  };

  for (std::size_t i = 0; i < m_cfg.maximumIterations; ++i) {
    auto maximumSharedMeasurements = *std::max_element(
        state.selectedTracks.begin(), state.selectedTracks.end(),
        sharedMeasurementsComperator);
    ACTS_VERBOSE(
        "maximum shared measurements "
        << state.sharedMeasurementsPerTrack[maximumSharedMeasurements]);
    if (state.sharedMeasurementsPerTrack[maximumSharedMeasurements] <
        m_cfg.maximumSharedHits) {
      break;
    }

    auto badTrack =
        *std::max_element(state.selectedTracks.begin(),
                          state.selectedTracks.end(), badTrackComperator);
    ACTS_VERBOSE("remove track " << badTrack);
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

  m_outputTrajectories(ctx, std::move(outputTrajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}
