// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/TrackFinding/detail/AmbiguityTrackClustering.hpp"

#include <map>
#include <unordered_map>
#include <vector>

#include "mlpack/methods/dbscan.hpp"

namespace Acts {

/// Clusterise tracks based on shared hits
///
/// @param trackMap : Multimap storing pair of track ID and vector of measurement ID. The keys are the number of measurement and are just there to focilitate the ordering.
/// @param tracks : Track container with all the track to be clustered
/// @param epsilon : Maximum distance between 2 tracks to be clustered
/// @param minPoints : Minimum number of tracks to create a cluster
/// @return an unordered map representing the clusters, the keys the ID of the primary track of each cluster and the store a vector of track IDs.
template <typename track_container_t, typename traj_t,
          template <typename> class holder_t>
std::unordered_map<int, std::vector<int>> dbscanTrackClustering(
    std::multimap<int, std::pair<int, std::vector<int>>>& trackMap,
    const Acts::TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    float epsilon = 0.07, int minPoints = 2) {
  // Unordered map associating a vector with all the track ID of a cluster to
  // the ID of the first track of the cluster
  std::unordered_map<int, std::vector<int>> cluster;
  // Unordered map associating hits to the ID of the first track of the
  // different clusters.
  std::unordered_map<int, int> hitToTrack;

  // DBSCAN algoritm from MLpack used in the track clustering
  mlpack::DBSCAN dbscan(epsilon, minPoints);

  arma::mat data(2, trackMap.size());
  int trackID = 0;
  arma::Row<size_t> assignments;

  // Get the input feature of the network for all the tracks
  for (const auto& [key, val] : trackMap) {
    auto traj = tracks.getTrack(val.first);
    data(0, trackID) = Acts::VectorHelpers::eta(traj.momentum());
    data(1, trackID) = Acts::VectorHelpers::phi(traj.momentum());
    trackID++;
  }
  size_t clusterNb = dbscan.Cluster(data, assignments);
  trackID = 0;

  // Cluster track with DBScan
  std::vector<std::multimap<int, std::pair<int, std::vector<int>>>>
      dbscanClusters(clusterNb);
  for (const auto& [key, val] : trackMap) {
    int clusterID = assignments(trackID);
    if (assignments(trackID) == SIZE_MAX) {
      cluster.emplace(val.first, std::vector<int>(1, val.first));
    } else {
      dbscanClusters[clusterID].emplace(key, val);
    }
    trackID++;
  }

  // Perform a subClustering of the DBScan cluster using the measurement ID
  // clustering
  for (const auto& dbscanCluster : dbscanClusters) {
    auto subCluster = Acts::detail::clusterDuplicateTracks(dbscanCluster);
    cluster.merge(subCluster);
    if (!subCluster.empty()) {
      std::cout << "Overlapping track ID, there must be an error" << std::endl;
    }
  }
  return cluster;
}

}  // namespace Acts
