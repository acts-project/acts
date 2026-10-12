// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/EventData/CalorimeterJet.hpp"
#include "ActsCalorimeter/Jets/CalorimeterJetReconstruction.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterMetadata.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"

#include <string>

#include <edm4hep/ClusterCollection.h>
#include <edm4hep/ReconstructedParticleCollection.h>

namespace ActsExamples {
/// Persist calorimeter jets as ReconstructedParticles with cluster relations.
/// PDG, charge and covariance remain unmeasured; the collection identifies
/// jets.
class EDM4hepCalorimeterJetOutputConverter final : public PodioOutputConverter {
 public:
  struct Config {
    std::string inputJets;
    std::string inputClusters;
    std::string inputHits;
    std::string inputEdmHits;
    std::string inputEdmClusters;
    std::string outputJets;
    /// Must match the preceding jet algorithm; validates constituent E-scheme
    /// sums.
    ActsCalorimeter::CalorimeterJetReconstruction::Config reconstruction;
  };
  explicit EDM4hepCalorimeterJetOutputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);
  ProcessCode execute(const AlgorithmContext& context) const override;
  std::vector<std::string> collections() const override;
  EDM4hepCalorimeterMetadata::Parameters metadata() const;
  const Config& config() const { return m_config; }

 private:
  Config m_config;
  ReadDataHandle<ActsCalorimeter::CalorimeterJetContainer> m_jets{this,
                                                                  "InputJets"};
  ReadDataHandle<ActsCalorimeter::CalorimeterClusterContainer> m_clusters{
      this, "InputClusters"};
  ReadDataHandle<ActsCalorimeter::CalorimeterHitContainer> m_hits{this,
                                                                  "InputHits"};
  PodioCollectionReadHandle<edm4hep::CalorimeterHitCollection> m_edmHits{
      this, "InputEdmHits"};
  PodioCollectionReadHandle<edm4hep::ClusterCollection> m_edmClusters{
      this, "InputEdmClusters"};
  PodioCollectionWriteHandle<edm4hep::ReconstructedParticleCollection> m_output{
      this, "OutputJets"};
};
}  // namespace ActsExamples
