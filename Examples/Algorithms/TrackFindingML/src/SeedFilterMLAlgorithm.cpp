// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingML/SeedFilterMLAlgorithm.hpp"

#include "Acts/Plugins/Mlpack/SeedFilterDBScanClustering.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iterator>
#include <map>

ActsExamples::SeedFilterMLAlgorithm::SeedFilterMLAlgorithm(
    ActsExamples::SeedFilterMLAlgorithm::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("SeedFilterMLAlgorithm", lvl),
      m_cfg(std::move(cfg)),
      m_seedClassifier(m_cfg.inputSeedFilterNN.c_str()) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.inputSimSeeds.empty()) {
    throw std::invalid_argument("Missing trajectories input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
  if (m_cfg.outputSimSeeds.empty()) {
    throw std::invalid_argument("Missing trajectories output collection");
  }
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputSimSeeds.initialize(m_cfg.inputSimSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputSimSeeds.initialize(m_cfg.outputSimSeeds);
}

ActsExamples::ProcessCode ActsExamples::SeedFilterMLAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Read input data
  const auto& seeds = m_inputSimSeeds(ctx);
  const auto& params = m_inputTrackParameters(ctx);
  if (seeds.size() != params.size()) {
    throw std::invalid_argument(
        "The number of seeds and parameters is different");
  }

  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      networkInput(seeds.size(), 14);
  std::vector<std::vector<double>> clusteringParams;
  std::vector<int> mapSeepIndex;
  // Loop over the seed and parameters to fill the input for the clustering
  // and the NN
  for (size_t i = 0; i < seeds.size(); i++) {
    double pT = std::abs(1.0 / params[i].parameters()[Acts::eBoundQOverP]) *
                std::sin(params[i].parameters()[Acts::eBoundTheta]);
    mapSeepIndex.push_back(i);
    size_t NNindex = mapSeepIndex.size() - 1;
    double eta =
        std::atanh(std::cos(params[i].parameters()[Acts::eBoundTheta]));
    double phi = params[i].parameters()[Acts::eBoundPhi];
    // Fill the clustering input
    clusteringParams.push_back({phi, eta, seeds[i].z() / 50, pT});
    // Fill the NN input
    networkInput.row(NNindex) << pT, eta, phi, seeds[i].sp()[0]->x(),
        seeds[i].sp()[0]->y(), seeds[i].sp()[0]->z(), seeds[i].sp()[1]->x(),
        seeds[i].sp()[1]->y(), seeds[i].sp()[1]->z(), seeds[i].sp()[2]->x(),
        seeds[i].sp()[2]->y(), seeds[i].sp()[2]->z(), seeds[i].z(),
        seeds[i].seedQuality();
  }

  // Cluster the tracks using DBscan
  auto cluster = Acts::dbscanSeedClustering(
      clusteringParams, m_cfg.epsilonDBScan, m_cfg.minPointsDBScan);

  // Select the ID of the track we want to keep
  std::vector<int> goodSeed =
      m_seedClassifier.solveAmbiguity(cluster, networkInput);

  // Create the output seed collection
  SimSeedContainer outputSeeds;
  outputSeeds.reserve(goodSeed.size());

  // Create the output track parameters collection
  TrackParametersContainer outputTrackParameters;
  outputTrackParameters.reserve(goodSeed.size());

  for (auto&& i : goodSeed) {
    outputSeeds.push_back(seeds[mapSeepIndex[i]]);
    outputTrackParameters.push_back(params[mapSeepIndex[i]]);
  }

  m_outputSimSeeds(ctx, SimSeedContainer{outputSeeds});
  m_outputTrackParameters(ctx, TrackParametersContainer{outputTrackParameters});

  return ActsExamples::ProcessCode::SUCCESS;
}
