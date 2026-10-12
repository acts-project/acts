// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/EventData/CalorimeterCluster.hpp"
#include "ActsCalorimeter/EventData/CalorimeterHit.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"

#include <cstdint>
#include <string>

#include <edm4hep/CalorimeterHitCollection.h>
#include <edm4hep/ClusterCollection.h>
#include <podio/UserDataCollection.h>

namespace ActsExamples {
/// Persist clusters and their hit relations, preserving native collection
/// order. Time and seed IDs use index-aligned sidecars; EDM4hep Cluster has no
/// such fields.
class EDM4hepCalorimeterClusterOutputConverter final
    : public PodioOutputConverter {
 public:
  struct Config {
    std::string inputClusters;
    std::string inputHits;
    std::string inputEdmHits;
    std::string outputClusters;
    std::string outputTimes;
    std::string outputSeedCellIds;
  };
  explicit EDM4hepCalorimeterClusterOutputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);
  ProcessCode execute(const AlgorithmContext& context) const override;
  std::vector<std::string> collections() const override;
  const Config& config() const { return m_config; }

 private:
  Config m_config;
  ReadDataHandle<ActsCalorimeter::CalorimeterClusterContainer> m_clusters{
      this, "InputClusters"};
  ReadDataHandle<ActsCalorimeter::CalorimeterHitContainer> m_hits{this,
                                                                  "InputHits"};
  PodioCollectionReadHandle<edm4hep::CalorimeterHitCollection> m_edmHits{
      this, "InputEdmHits"};
  PodioCollectionWriteHandle<edm4hep::ClusterCollection> m_output{
      this, "OutputClusters"};
  PodioCollectionWriteHandle<podio::UserDataCollection<double>> m_times{
      this, "OutputTimes"};
  PodioCollectionWriteHandle<podio::UserDataCollection<std::uint64_t>> m_seeds{
      this, "OutputSeeds"};
};
}  // namespace ActsExamples
