// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Calorimeter/EDM4hepCalorimeterMetadata.hpp"

#include "ExampleFixture.hpp"

namespace ActsExamples::CalorimeterFixture {
inline EDM4hepCalorimeterInputConverter::Config inputConfig() {
  return {"events", "SimCaloHits", "calo_deposits", "calo_sources", cellCentre};
}
inline CalorimeterDigitizationAlgorithm::Config responseConfig() {
  CalorimeterDigitizationAlgorithm::Config config;
  config.inputSimHits = "calo_deposits";
  config.outputHits = "calo_hits";
  config.response.energyScale = 2;
  config.response.energyThreshold = 0.01 * Acts::UnitConstants::GeV;
  config.response.timeMin = 0;
  config.response.timeMax = 10 * Acts::UnitConstants::ns;
  return config;
}
inline CalorimeterClusteringAlgorithm::Config clusteringConfig() {
  CalorimeterClusteringAlgorithm::Config config;
  config.inputHits = "calo_hits";
  config.outputClusters = "calo_clusters";
  config.clustering.seedEnergyThreshold = 0.03 * Acts::UnitConstants::GeV;
  config.clustering.neighbourEnergyThreshold = 0.01 * Acts::UnitConstants::GeV;
  config.neighbours = {{highCellId, highCellId + 1}};
  return config;
}
inline EDM4hepCalorimeterOutputConverter::Config hitOutputConfig() {
  return {"calo_hits", "calo_deposits", "calo_sources", "CaloHits",
          "CaloLinks"};
}
inline EDM4hepCalorimeterClusterOutputConverter::Config clusterOutputConfig() {
  return {"calo_clusters", "calo_hits",        "CaloHits",
          "CaloClusters",  "CaloClusterTimes", "CaloClusterSeedCellIds"};
}
inline EDM4hepCalorimeterMetadata::Config metadataConfig() {
  EDM4hepCalorimeterMetadata::Config config;
  config.outputFrame = "calo_events";
  config.geometryIdentifier =
      "synthetic calorimeter IO fixture v2; global cell centres in mm";
  config.input = inputConfig();
  config.response = responseConfig();
  config.clustering = clusteringConfig();
  config.hitOutput = hitOutputConfig();
  config.clusterOutput = clusterOutputConfig();
  return config;
}
}  // namespace ActsExamples::CalorimeterFixture
