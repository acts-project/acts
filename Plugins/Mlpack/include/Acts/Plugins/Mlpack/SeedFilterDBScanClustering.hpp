// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <unordered_map>
#include <vector>

#include "mlpack/methods/dbscan.hpp"

namespace Acts {

/// Clusterise seed based on their Z position, their direction and their
/// momentum using DBScan
///
/// @param input : Input parameters for the clustering (phi, eta, z, Pt/10)
/// @param epsilon : Maximum distance between 2 tracks to be clustered
/// @param minPoints : Minimum number of tracks to create a cluster
/// @return an unordered map representing the clusters, the keys the ID of the primary track of each cluster and the store a vector of track IDs.
std::vector<std::vector<int>> dbscanSeedClustering(
    const std::vector<std::vector<double>>& input, float epsilon = 0.07,
    int minPoints = 2) {
  // DBSCAN algoritm from MLpack used in the seed clustering
  mlpack::DBSCAN dbscan(epsilon, minPoints);

  // Compute the space dimension of the input
  int dim = input[0].size();

  arma::mat data(dim, input.size());
  arma::Row<size_t> assignments;

  size_t trackID = 0;
  // Get the input feature of the network for all the tracks
  for (const auto& param : input) {
    for (int i = 0; i < dim; i++) {
      data(i, trackID) = param[i];
    }
    trackID++;
  }
  size_t clusterNb = dbscan.Cluster(data, assignments);

  // Cluster track with DBScan
  std::vector<std::vector<int>> cluster(clusterNb);
  for (size_t iD = 0; iD < input.size(); iD++) {
    int clusterID = assignments(iD);
    if (assignments(iD) == SIZE_MAX) {
      cluster.push_back(std::vector<int>(1, iD));
    } else {
      cluster[clusterID].push_back(iD);
    }
  }
  return cluster;
}

}  // namespace Acts
