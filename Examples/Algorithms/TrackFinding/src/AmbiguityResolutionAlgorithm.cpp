// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/AmbiguityResolutionAlgorithm.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iterator>
#include <numeric>
#include <stdexcept>

ActsExamples::AmbiguityResolutionAlgorithm::AmbiguityResolutionAlgorithm(
    ActsExamples::AmbiguityResolutionAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("AmbiguityResolutionAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links input collection");
  }
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
}

namespace {

// TODO this is somewhat duplicated in TrackFindingAlgorithm.hpp
// TODO we should make a common implementation in the core at some point
std::vector<std::size_t> computeSharedHits(
    const ActsExamples::IndexSourceLinkContainer& sourceLinks,
    const ActsExamples::TrajectoriesContainer& trajectories,
    const std::vector<uint32_t>& trackIndices,
    const std::vector<std::pair<size_t, size_t>>& trackTips) {
  std::vector<std::size_t> hitCountPerMeasurement(sourceLinks.size(), 0);

  for (auto indexTrack : trackIndices) {
    const auto [indexTraj, tip] = trackTips[indexTrack];
    const auto& traj = trajectories[indexTraj];

    traj.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
      if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }

      const std::size_t indexHit =
          state.uncalibratedSourceLink()
              .template get<ActsExamples::IndexSourceLink>()
              .index();

      ++hitCountPerMeasurement[indexHit];

      return true;
    });
  }

  std::vector<std::size_t> sharedHitCountPerTrack(trackIndices.size(), 0);

  for (std::size_t i = 0; i < trackIndices.size(); ++i) {
    const auto indexTrack = trackIndices[i];
    const auto [indexTraj, tip] = trackTips[indexTrack];
    const auto& traj = trajectories[indexTraj];

    traj.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
      if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }

      const std::size_t indexHit =
          state.uncalibratedSourceLink()
              .template get<ActsExamples::IndexSourceLink>()
              .index();

      if (hitCountPerMeasurement[indexHit] > 1) {
        ++sharedHitCountPerTrack[i];
      }

      return true;
    });
  }

  return sharedHitCountPerTrack;
}

std::size_t computeTrackHits(const Acts::VectorMultiTrajectory& multiTrajectory,
                             const std::size_t tip) {
  std::size_t result = 0;

  multiTrajectory.visitBackwards(tip, [&](const auto&) { ++result; });

  return result;
}

}  // namespace

ActsExamples::ProcessCode ActsExamples::AmbiguityResolutionAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  TrackParametersContainer trackParameters;
  std::vector<std::pair<size_t, size_t>> trackTips;

  for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
    const auto& traj = trajectories[iTraj];
    for (auto tip : traj.tips()) {
      if (!traj.hasTrackParameters(tip)) {
        continue;
      }
      trackParameters.push_back(traj.trackParameters(tip));
      trackTips.emplace_back(iTraj, tip);
    }
  }

  std::vector<uint32_t> hitCount(trackParameters.size(), 0);
  for (std::size_t i = 0; i < trackParameters.size(); ++i) {
    const auto [iTraj, tip] = trackTips[i];
    const auto& traj = trajectories[iTraj];
    hitCount[i] = computeTrackHits(traj.multiTrajectory(), tip);
  }

  std::vector<uint32_t> trackIndices(trackParameters.size());
  std::iota(std::begin(trackIndices), std::end(trackIndices), 0);

  while (true) {
    const auto sharedHits =
        computeSharedHits(sourceLinks, trajectories, trackIndices, trackTips);

    if (sharedHits.empty() ||
        *std::max_element(std::begin(sharedHits), std::end(sharedHits)) <
            m_cfg.maximumSharedHits) {
      break;
    }

    std::vector<float> relativeSharedHits(trackIndices.size(), 0);
    for (std::size_t i = 0; i < trackIndices.size(); ++i) {
      const auto indexTrack = trackIndices[i];
      relativeSharedHits[i] = 1.0f * sharedHits[i] / hitCount[indexTrack];
    }

    const auto maxRelativeSharedHits = std::max_element(
        std::begin(relativeSharedHits), std::end(relativeSharedHits));
    const auto index =
        std::distance(std::begin(relativeSharedHits), maxRelativeSharedHits);
    trackIndices.erase(std::begin(trackIndices) + index);
  }

  if (trackIndices.size() == trackParameters.size()) {
    const auto sharedHits =
        computeSharedHits(sourceLinks, trajectories, trackIndices, trackTips);

    std::vector<float> relativeSharedHits(trackIndices.size(), 0);
    for (std::size_t i = 0; i < trackIndices.size(); ++i) {
      const auto indexTrack = trackIndices[i];
      relativeSharedHits[i] = 1.0f * sharedHits[i] / hitCount[indexTrack];
    }
  }

  ACTS_INFO("Resolved to " << trackIndices.size() << " tracks from "
                           << trackParameters.size());

  TrajectoriesContainer outputTrajectories;
  outputTrajectories.reserve(trajectories.size());
  for (std::size_t iTraj = 0; iTraj < trajectories.size(); ++iTraj) {
    const auto& traj = trajectories[iTraj];

    std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;
    Trajectories::IndexedParameters parameters;

    for (auto iTrack : trackIndices) {
      if (trackTips[iTrack].first != iTraj) {
        continue;
      }
      const auto tip = trackTips[iTrack].second;
      tips.push_back(tip);
      parameters.emplace(tip, trackParameters[iTrack]);
    }

    outputTrajectories.emplace_back(traj.multiTrajectoryPtr(), tips,
                                    parameters);
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(outputTrajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}
