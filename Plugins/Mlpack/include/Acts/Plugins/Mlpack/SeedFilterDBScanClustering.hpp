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

/// Clusters seed based on their direction, their Z impact parameter and their
/// momentum using DBScan
///
/// @param input : Input parameters for the clustering (phi, eta, z, Pt)
/// @param epsilon : Maximum distance between 2 tracks to be clustered
/// @param minPoints : Minimum number of tracks to create a cluster
/// @return an unordered map representing the clusters, the keys the ID of the primary seed of each cluster and the stored value a vector of seed IDs.
std::vector<std::vector<std::size_t>> dbscanSeedClustering(
    const std::vector<std::vector<double>>& input, float epsilon = 0.03,
    int minPoints = 2) {
  // DBSCAN algorithm from MLpack used in the seed clustering
  mlpack::DBSCAN dbscan(epsilon, minPoints);

  // Compute the space dimension of the input
  int dim = input[0].size();

  // Prepare the input for the DBScan
  arma::mat data(dim, input.size());
  arma::Row<std::size_t> assignments;
  std::size_t trackID = 0;
  for (const auto& param : input) {
    for (int i = 0; i < dim; i++) {
      data(i, trackID) = param[i];
    }
    trackID++;
  }
  // Cluster track with DBScan
  std::size_t clusterNb = dbscan.Cluster(data, assignments);

  // Prepare the output
  std::vector<std::vector<std::size_t>> cluster(clusterNb,
                                                std::vector<std::size_t>());
  for (std::size_t iD = 0; iD < input.size(); iD++) {
    std::size_t clusterID = assignments(iD);
    if (assignments(iD) == SIZE_MAX) {
      cluster.push_back(std::vector<std::size_t>(1, iD));
    } else {
      cluster[clusterID].push_back(iD);
    }
  }
  return cluster;
}

}  // namespace Acts
