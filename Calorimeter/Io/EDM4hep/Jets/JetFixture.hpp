// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Calorimeter/CalorimeterJetAlgorithm.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterJetOutputConverter.hpp"

#include "../ReconstructionFixture.hpp"

namespace ActsExamples::CalorimeterFixture {
inline CalorimeterJetAlgorithm::Config jetConfig() {
  CalorimeterJetAlgorithm::Config config;
  config.inputClusters = "calo_clusters";
  config.outputJets = "calo_jets";
  return config;
}
inline EDM4hepCalorimeterJetOutputConverter::Config jetOutputConfig() {
  const auto jets = jetConfig();
  return {jets.outputJets, jets.inputClusters, "calo_hits", "CaloHits",
          "CaloClusters",  "CaloJets",         jets.jets};
}
}  // namespace ActsExamples::CalorimeterFixture
