// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsCalorimeter/EventData/CalorimeterHit.hpp"
#include "ActsExamples/Calorimeter/EDM4hepCalorimeterSource.hpp"
#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioOutputConverter.hpp"

#include <string>

#include <edm4hep/CaloHitSimCaloHitLinkCollection.h>
#include <edm4hep/CalorimeterHitCollection.h>

namespace ActsExamples {

/// Write calibrated hits and energy-fraction links to original simulated hits.
class EDM4hepCalorimeterOutputConverter final : public PodioOutputConverter {
 public:
  struct Config {
    std::string inputHits;
    std::string inputDeposits;
    std::string inputSources;
    std::string outputHits;
    std::string outputLinks;
  };

  explicit EDM4hepCalorimeterOutputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  ProcessCode execute(const AlgorithmContext& context) const override;
  std::vector<std::string> collections() const override;
  const Config& config() const { return m_config; }

 private:
  Config m_config;
  ReadDataHandle<ActsCalorimeter::CalorimeterHitContainer> m_hits{this,
                                                                  "InputHits"};
  ReadDataHandle<ActsCalorimeter::SimCalorimeterHitContainer> m_deposits{
      this, "InputDeposits"};
  ReadDataHandle<EDM4hepCalorimeterSourceContainer> m_sources{this,
                                                              "InputSources"};
  PodioCollectionWriteHandle<edm4hep::CalorimeterHitCollection> m_outputHits{
      this, "OutputHits"};
  PodioCollectionWriteHandle<edm4hep::CaloHitSimCaloHitLinkCollection>
      m_outputLinks{this, "OutputLinks"};
};

}  // namespace ActsExamples
