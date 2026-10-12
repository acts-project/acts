// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/EDM4hepCalorimeterInputConverter.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <edm4hep/SimCalorimeterHitCollection.h>
#include <podio/Frame.h>

namespace ActsExamples {

EDM4hepCalorimeterInputConverter::EDM4hepCalorimeterInputConverter(
    const Config& config, Acts::Logging::Level level)
    : PodioInputConverter(
          "EDM4hepCalorimeterInput", config.inputFrame,
          Acts::getDefaultLogger("EDM4hepCalorimeterInput", level)),
      m_config(config) {
  if (config.inputFrame.empty() || config.inputSimHits.empty() ||
      config.outputDeposits.empty() || config.outputSources.empty() ||
      config.outputDeposits == config.outputSources ||
      config.inputFrame == config.outputDeposits ||
      config.inputFrame == config.outputSources || !config.cellCentre) {
    throw std::invalid_argument(
        "Invalid calorimeter input converter configuration");
  }
  m_deposits.initialize(config.outputDeposits);
  m_sources.initialize(config.outputSources);
}

ProcessCode EDM4hepCalorimeterInputConverter::convert(
    const AlgorithmContext& context, const podio::Frame& frame) const {
  const auto& input =
      frame.get<edm4hep::SimCalorimeterHitCollection>(m_config.inputSimHits);
  ActsCalorimeter::SimCalorimeterHitContainer deposits;
  EDM4hepCalorimeterSourceContainer sources;
  std::map<std::uint64_t, Acts::Vector3> centres;
  std::set<std::pair<std::uint32_t, int>> usedContributions;
  for (const auto& hit : input) {
    if (!hit.isAvailable() || !std::isfinite(hit.getEnergy()) ||
        hit.getEnergy() < 0) {
      throw std::invalid_argument("Invalid EDM4hep simulated calorimeter hit");
    }
    const auto id = hit.getCellID();
    auto [centre, inserted] = centres.try_emplace(id, Acts::Vector3::Zero());
    if (inserted) {
      centre->second = m_config.cellCentre(id);
      if (!centre->second.allFinite()) {
        throw std::invalid_argument("Nonfinite calorimeter cell centre");
      }
    }
    double energySum = 0;
    for (const auto& contribution : hit.getContributions()) {
      if (!contribution.isAvailable() ||
          !std::isfinite(contribution.getEnergy()) ||
          contribution.getEnergy() < 0 ||
          !std::isfinite(contribution.getTime())) {
        throw std::invalid_argument("Invalid EDM4hep calorimeter contribution");
      }
      const auto objectId = contribution.getObjectID();
      if (objectId.index < 0 ||
          !usedContributions.emplace(objectId.collectionID, objectId.index)
               .second) {
        throw std::invalid_argument(
            "Repeated or unregistered calorimeter contribution");
      }
      energySum += contribution.getEnergy();
      deposits.push_back({id, centre->second,
                          contribution.getEnergy() * Acts::UnitConstants::GeV,
                          contribution.getTime() * Acts::UnitConstants::ns});
      sources.push_back({hit, contribution});
    }
    // A positive hit without contributions has neither deposit times nor
    // traceable deposits. Do not invent a zero time or rescale contributions.
    const double tolerance =
        16 * std::numeric_limits<float>::epsilon() *
        std::max(energySum, static_cast<double>(hit.getEnergy()));
    if (!std::isfinite(energySum) ||
        std::abs(energySum - hit.getEnergy()) > tolerance) {
      throw std::invalid_argument(
          "Simulated hit energy disagrees with its contributions");
    }
  }
  m_deposits(context, std::move(deposits));
  m_sources(context, std::move(sources));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
