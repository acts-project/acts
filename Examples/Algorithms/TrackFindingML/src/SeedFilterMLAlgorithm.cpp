// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFindingML/SeedFilterMLAlgorithm.hpp"

#include "ActsExamples/TrackFindingML/SeedFilterDBScanClustering.hpp"

using namespace Acts;
using namespace ActsPlugins;

namespace ActsExamples {

SeedFilterMLAlgorithm::SeedFilterMLAlgorithm(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("SeedFilterMLAlgorithm", std::move(logger)),
      m_cfg(cfg),
      m_seedClassifier(m_cfg.inputSeedFilterNN.c_str()) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters input collection");
  }
  if (m_cfg.inputSeeds.empty()) {
    throw std::invalid_argument("Missing seed input collection");
  }
  if (m_cfg.outputTrackParameters.empty()) {
    throw std::invalid_argument("Missing track parameters output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seed output collection");
  }
  m_inputTrackParameters.initialize(m_cfg.inputTrackParameters);
  m_inputSeeds.initialize(m_cfg.inputSeeds);
  m_outputTrackParameters.initialize(m_cfg.outputTrackParameters);
  m_outputSeeds.initialize(m_cfg.outputSeeds);
}

ProcessCode SeedFilterMLAlgorithm::execute(const AlgorithmContext& ctx) const {
  // Read input data
  const auto& seeds = m_inputSeeds(ctx);
  const auto& params = m_inputTrackParameters(ctx);
  if (seeds.size() != params.size()) {
    throw std::invalid_argument(
        "The number of seeds and track parameters is different");
  }

  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      networkInput(seeds.size(), 14);
  std::vector<std::array<double, 4>> clusteringParams;
  // Loop over the seed and parameters to fill the input for the clustering
  // and the NN
  for (std::size_t i = 0; i < seeds.size(); i++) {
    // Compute the track parameters
    double pT = std::abs(1.0 / params[i].parameters()[eBoundQOverP]) *
                std::sin(params[i].parameters()[eBoundTheta]);
    double eta = std::atanh(std::cos(params[i].parameters()[eBoundTheta]));
    double phi = params[i].parameters()[eBoundPhi];

    // Fill and weight the clustering inputs
    clusteringParams.push_back({phi / m_cfg.clusteringWeighPhi,
                                eta / m_cfg.clusteringWeighEta,
                                seeds[i].vertexZ() / m_cfg.clusteringWeighZ,
                                pT / m_cfg.clusteringWeighPt});

    // Fill the NN input
    networkInput.row(i) << pT, eta, phi, seeds[i].spacePoints()[0].x(),
        seeds[i].spacePoints()[0].y(), seeds[i].spacePoints()[0].z(),
        seeds[i].spacePoints()[1].x(), seeds[i].spacePoints()[1].y(),
        seeds[i].spacePoints()[1].z(), seeds[i].spacePoints()[2].x(),
        seeds[i].spacePoints()[2].y(), seeds[i].spacePoints()[2].z(),
        seeds[i].vertexZ(), seeds[i].quality();
  }

  // Cluster the tracks using DBscan
  auto cluster = dbscanSeedClustering(clusteringParams, m_cfg.epsilonDBScan,
                                      m_cfg.minPointsDBScan);

  // Select the ID of the track we want to keep
  std::vector<std::size_t> goodSeedIndices = m_seedClassifier.solveAmbiguity(
      cluster, networkInput, m_cfg.minSeedScore);

  // Create the output seed collection
  SeedContainer outputSeeds;
  outputSeeds.assignSpacePointContainer(seeds.spacePointContainer());
  outputSeeds.reserve(goodSeedIndices.size());

  // Create the output track parameters collection
  TrackParametersContainer outputTrackParameters;
  outputTrackParameters.reserve(goodSeedIndices.size());

  for (std::size_t i : goodSeedIndices) {
    auto newSeed = outputSeeds.createSeed();
    newSeed.assignSpacePointIndices(seeds[i].spacePointIndices());
    newSeed.vertexZ() = seeds[i].vertexZ();
    newSeed.quality() = seeds[i].quality();
    outputTrackParameters.push_back(params[i]);
  }

  m_outputSeeds(ctx, std::move(outputSeeds));
  m_outputTrackParameters(ctx, std::move(outputTrackParameters));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
