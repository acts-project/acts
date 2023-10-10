// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/detail/AmbiguityTrackClustering.hpp"

#include <iterator>

std::unordered_map<int, std::vector<int>> Acts::detail::clusterDuplicateTracks(
    const std::multimap<int, std::pair<int, std::vector<int>>>& trackMap) {
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
      cluster.emplace(track->second.first,
                      std::vector<int>(1, track->second.first));
      for (const auto& hit : hits) {
        // Add the hits of the new cluster to the hitToTrack
        hitToTrack.emplace(hit, track->second.first);
      }
    }
  }
  return cluster;
}
