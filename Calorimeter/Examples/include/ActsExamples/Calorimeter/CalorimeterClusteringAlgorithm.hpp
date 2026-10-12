// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/Clustering/CalorimeterClusterer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <string>

namespace ActsExamples {

/// ActsExamples adapter with fixed geometry adjacency and event-owned clusters.
class CalorimeterClusteringAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputHits;
    std::string outputClusters;
    ActsCalorimeter::CalorimeterClusterer::Config clustering;
    /// Snapshot of geometry adjacency, copied at construction. This is not an
    /// event collection and no topology is inferred from the input hits.
    ActsCalorimeter::CellNeighbourContainer neighbours;
  };

  explicit CalorimeterClusteringAlgorithm(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  ProcessCode execute(const AlgorithmContext& context) const override;
  const Config& config() const { return m_config; }

 private:
  Config m_config;
  ActsCalorimeter::CalorimeterClusterer m_clusterer;
  ReadDataHandle<ActsCalorimeter::CalorimeterHitContainer> m_input{this,
                                                                   "InputHits"};
  WriteDataHandle<ActsCalorimeter::CalorimeterClusterContainer> m_output{
      this, "OutputClusters"};
};

}  // namespace ActsExamples
