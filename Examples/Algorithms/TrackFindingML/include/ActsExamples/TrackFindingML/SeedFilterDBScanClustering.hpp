// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/DBScan.hpp"

#include <map>
#include <unordered_map>
#include <vector>

namespace Acts {

/// Clusters seed based on their direction, their Z impact parameter and their
/// momentum using DBScan
///
/// @param input : Input parameters for the clustering (phi, eta, z, Pt)
/// @param epsilon : Maximum distance between 2 tracks to be clustered
/// @param minPoints : Minimum number of seed to create a cluster
/// @return an unordered map representing the clusters, the keys the ID of the primary seed of each cluster and the stored value a vector of seed IDs.
std::vector<std::vector<std::size_t>> dbscanSeedClustering(
    const std::vector<std::array<double, 4>>& input, float epsilon = 0.03,
    int minPoints = 2) {
  // Initialize a DBScan of dimension 4 (phi, eta, z, Pt)
  using DBSCAN = Acts::DBScan<4, double, 4>;
  DBSCAN dbscan(epsilon, minPoints, true);

  // create a chrono to measure the time
  std::vector<int> clusterAssignments;
  std::size_t clusterNb = dbscan.cluster(input, clusterAssignments);

  // Prepare the output
  std::vector<std::vector<std::size_t>> cluster(clusterNb,
                                                std::vector<std::size_t>());
  for (std::size_t iD = 0; iD < input.size(); iD++) {
    int clusterID = clusterAssignments[iD];
    cluster[clusterID].push_back(iD);
  }

  return cluster;
}

}  // namespace Acts
