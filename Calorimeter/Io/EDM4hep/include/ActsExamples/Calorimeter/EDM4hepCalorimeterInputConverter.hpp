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
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/Podio/PodioInputConverter.hpp"

#include <cstdint>
#include <functional>
#include <string>

namespace ActsExamples {

/// Expand EDM4hep calorimeter contributions into native deposits in ACTS units.
class EDM4hepCalorimeterInputConverter final : public PodioInputConverter {
 public:
  struct Config {
    std::string inputFrame = "events";
    std::string inputSimHits;
    std::string outputDeposits;
    std::string outputSources;
    /// Required geometry lookup returning a global cell centre in ACTS units.
    /// EDM4hep hit/step positions are not assumed to be cell centres. Called
    /// once per cell per event; must support concurrent calls for parallel
    /// runs.
    std::function<Acts::Vector3(std::uint64_t)> cellCentre;
  };

  explicit EDM4hepCalorimeterInputConverter(
      const Config& config, Acts::Logging::Level level = Acts::Logging::INFO);

  const Config& config() const { return m_config; }

 private:
  ProcessCode convert(const AlgorithmContext& context,
                      const podio::Frame& frame) const override;

  Config m_config;
  WriteDataHandle<ActsCalorimeter::SimCalorimeterHitContainer> m_deposits{
      this, "OutputDeposits"};
  WriteDataHandle<EDM4hepCalorimeterSourceContainer> m_sources{this,
                                                               "OutputSources"};
};

}  // namespace ActsExamples
