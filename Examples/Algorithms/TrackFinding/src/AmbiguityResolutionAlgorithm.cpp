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
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing track parameter tips input collection");
  }
  if (m_cfg.inputTrackParametersTips.empty()) {
    throw std::invalid_argument("Missing track parameters input collection");
  }
  if (m_cfg.outputTrackIndices.empty()) {
    throw std::invalid_argument("Missing track indices output collection");
  }
}

namespace {

std::vector<std::size_t> computeSharedHits(
    const ActsExamples::IndexSourceLinkContainer& sourceLinks,
    const ActsExamples::TrajectoriesContainer& trajectories,
    const std::vector<uint32_t>& trackIndices,
    const std::vector<std::pair<size_t, size_t>>& trackTips) {
  std::vector<std::size_t> hitCountPerMeasurement(sourceLinks.size(), 0);

  for (std::size_t i = 0; i < trackIndices.size(); ++i) {
    const auto indexTrack = trackIndices[i];
    const auto [indexTraj, tip] = trackTips[indexTrack];
    const auto& traj = trajectories[indexTraj];

    traj.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
      if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }

      const std::size_t indexHit =
          static_cast<const ActsExamples::IndexSourceLink&>(
              state.uncalibrated())
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
          static_cast<const ActsExamples::IndexSourceLink&>(
              state.uncalibrated())
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
  const auto& trackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  const auto& trackTips =
      ctx.eventStore.get<std::vector<std::pair<size_t, size_t>>>(
          m_cfg.inputTrackParametersTips);

  std::vector<uint32_t> hitCount(trackParameters.size(), 0);
  for (std::size_t i = 0; i < trackParameters.size(); ++i) {
    const auto [indexTraj, tip] = trackTips[i];
    const auto& traj = trajectories[indexTraj];

    hitCount[i] = computeTrackHits(traj.multiTrajectory(), tip);
  }

  std::vector<uint32_t> trackIndices(trackParameters.size());
  std::iota(std::begin(trackIndices), std::end(trackIndices), 0);

  while (true) {
    const auto sharedHits =
        computeSharedHits(sourceLinks, trajectories, trackIndices, trackTips);

    std::vector<float> relativeSharedHits(trackIndices.size(), 0);
    for (std::size_t i = 0; i < trackIndices.size(); ++i) {
      const auto indexTrack = trackIndices[i];

      relativeSharedHits[i] = 1.0f * sharedHits[i] / hitCount[indexTrack];
    }

    auto it = std::max_element(std::begin(relativeSharedHits),
                               std::end(relativeSharedHits));
    if (it == std::end(relativeSharedHits)) {
      break;
    }

    const auto index = it - std::begin(relativeSharedHits);
    if (sharedHits[index] < this->m_cfg.maximumSharedHits) {
      break;
    }

    trackIndices.erase(std::begin(trackIndices) + index);
  }

  ACTS_INFO("Resolved to " << trackIndices.size() << " tracks from "
                           << trackParameters.size());

  TrackParametersContainer outputTrackParameters;
  outputTrackParameters.reserve(trackIndices.size());
  for (std::size_t i = 0; i < trackIndices.size(); ++i) {
    outputTrackParameters.push_back(trackParameters[trackIndices[i]]);
  }

  ctx.eventStore.add(m_cfg.outputTrackParameters,
                     std::move(outputTrackParameters));
  ctx.eventStore.add(m_cfg.outputTrackIndices, std::move(trackIndices));
  return ActsExamples::ProcessCode::SUCCESS;
}
