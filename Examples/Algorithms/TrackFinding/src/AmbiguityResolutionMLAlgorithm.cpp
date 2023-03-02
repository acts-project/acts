// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/AmbiguityResolutionMLAlgorithm.hpp"

#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <chrono>
#include <iterator>
#include <map>
#include <numeric>
#include <stdexcept>

#include <Eigen/Dense>
#include <core/session/onnxruntime_cxx_api.h>

ActsExamples::AmbiguityResolutionMLAlgorithm::AmbiguityResolutionMLAlgorithm(
    ActsExamples::AmbiguityResolutionMLAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("AmbiguityResolutionMLAlgorithm", lvl),
      m_cfg(std::move(cfg)),
      m_env(ORT_LOGGING_LEVEL_WARNING, "MLClassifier"),
      m_duplicateClassifier(m_env, m_cfg.inputDuplicateNN.c_str()) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
}

namespace {

/// Clusterise tracks based on shared hits
///
/// @param trackMap : Multimap storing pair of track ID and vector of measurement ID. The keys are the number of measurement and are just there to focilitate the ordering.
/// @return an unordered map representing the clusters, the keys the ID of the primary track of each cluster and the store a vector of track IDs.
std::unordered_map<int, std::vector<int>> clusterTracks(
    std::multimap<int, std::pair<int, std::vector<int>>> trackMap) {
  // Unordered map associating a vector with all the track ID of a cluster to
  // the ID of the first track of the cluster
  std::unordered_map<int, std::vector<int>> cluster;
  // Unordered map associating hits to the ID of the first track of the
  // different clusters.
  std::unordered_map<int, int> hitToTrack;

  // Loop over all the tracks
  for (auto track = trackMap.rbegin(); track != trackMap.rend(); ++track) {
    std::vector<int> hits = track->second.second;
    auto matchedTrack = hitToTrack.end();
    // Loop over all the hits in the track
    for (auto hit = hits.begin(); hit != hits.end(); hit++) {
      // Check if the hit is already associated to a track
      matchedTrack = hitToTrack.find(*hit);
      if (matchedTrack != hitToTrack.end()) {
        // Add the track to the cluster associated to the matched track
        cluster.at(matchedTrack->second).push_back(track->second.first);
        break;
      }
    }
    // None of the hits have been matched to a track create a new cluster
    if (matchedTrack == hitToTrack.end()) {
      cluster.emplace(track->second.first, {track->second.first});
      for (const auto& hit : hits) {
        // Add the hits of the new cluster to the hitToTrack
        hitToTrack.emplace(hit, track->second.first);
      }
    }
  }
  return cluster;
}
}  // namespace

ActsExamples::ProcessCode ActsExamples::AmbiguityResolutionMLAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& trajectories =
      ctx.eventStore.get<TrajectoriesContainer>(m_cfg.inputTrajectories);

  TrackParametersContainer trackParameters;
  std::vector<std::pair<int, int>> trackTips;
  int trackID = 0;
  int iTraj = 0;
  std::multimap<int, std::pair<int, std::vector<int>>> trackMap;

  // Loop over all the trajectories in the events
  for (const auto& traj : trajectories) {
    for (auto tip : traj.tips()) {
      if (!traj.hasTrackParameters(tip)) {
        continue;
      }
      // Store the hits id for the trajectory and compute the number of
      // measurement
      std::vector<int> hits;
      int nbMeasurements = 0;
      traj.multiTrajectory().visitBackwards(tip, [&](const auto& state) {
        if (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          int indexHit = state.uncalibratedSourceLink()
                             .template get<ActsExamples::IndexSourceLink>()
                             .index();
          hits.emplace_back(indexHit);
          ++nbMeasurements;
        }
      });
      if (nbMeasurements < m_cfg.nMeasurementsMin) {
        continue;
      }
      trackMap.emplace(nbMeasurements, std::make_pair(trackID, hits));
      auto param = traj.trackParameters(tip);
      trackParameters.emplace_back(param);
      trackTips.emplace_back(iTraj, tip);
      trackID++;
    }
    iTraj++;
  }
  // Performe the share hit based clustering
  auto clusters = clusterTracks(trackMap);
  Acts::NetworkBatchInput networkInput(trackID + 1, 8);
  trackID = 0;
  // Get the input feature of the network for all the tracks
  for (const auto& [key, val] : clusters) {
    for (const auto& track : val) {
      std::pair<int, int> tips = trackTips.at(track);
      TrackParameters parameters = trackParameters.at(track);
      auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(
          trajectories[tips.first].multiTrajectory(), tips.second);
      networkInput(trackID, 0) = trajState.nStates;
      networkInput(trackID, 1) = trajState.nMeasurements;
      networkInput(trackID, 2) = trajState.nOutliers;
      networkInput(trackID, 3) = trajState.nHoles;
      networkInput(trackID, 4) = trajState.NDF;
      networkInput(trackID, 5) = (trajState.chi2Sum * 1.0) / trajState.NDF;
      networkInput(trackID, 6) =
          Acts::VectorHelpers::eta(parameters.momentum());
      networkInput(trackID, 7) =
          Acts::VectorHelpers::phi(parameters.momentum());
      trackID++;
    }
  }
  // Use the network to compute a score for all the tracks.
  std::vector<std::vector<float>> outputTensor =
      m_duplicateClassifier.runONNXInference(networkInput);
  std::vector<int> goodTracks;
  int iOut = 0;

  // Loop over all the cluster and only keep the track with the highest score in
  // each cluster
  for (const auto& [key, val] : clusters) {
    int bestTrackID = 0;
    float bestTrackScore = 0;
    for (const auto& track : val) {
      if (outputTensor[iOut][0] > bestTrackScore) {
        bestTrackScore = outputTensor[iOut][0];
        bestTrackID = track;
      }
      iOut++;
    }
    goodTracks.push_back(bestTrackID);
  }

  // Create an output multitrajectory based of the good tracks
  TrajectoriesContainer outputTrajectories;
  outputTrajectories.reserve(goodTracks.size());
  for (auto&& iTrack : goodTracks) {
    const auto& outputTips = trackTips.at(iTrack);

    std::vector<Acts::MultiTrajectoryTraits::IndexType> tips;
    Trajectories::IndexedParameters parameters;
    tips.push_back(outputTips.second);
    parameters.emplace(outputTips.second, trackParameters[iTrack]);

    outputTrajectories.emplace_back(
        trajectories[outputTips.first].multiTrajectory(), tips, parameters);
  }

  // Add our output multitrajectories to the event store
  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(outputTrajectories));

  return ActsExamples::ProcessCode::SUCCESS;
}
