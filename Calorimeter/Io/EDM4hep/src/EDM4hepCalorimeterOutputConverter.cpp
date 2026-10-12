// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Calorimeter/EDM4hepCalorimeterOutputConverter.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace ActsExamples {
namespace {

float checkedFloat(double value) {
  if (!std::isfinite(value) ||
      std::abs(value) > std::numeric_limits<float>::max()) {
    throw std::overflow_error("Calorimeter value exceeds EDM4hep float range");
  }
  return static_cast<float>(value);
}

}  // namespace

EDM4hepCalorimeterOutputConverter::EDM4hepCalorimeterOutputConverter(
    const Config& config, Acts::Logging::Level level)
    : PodioOutputConverter(
          "EDM4hepCalorimeterOutput",
          Acts::getDefaultLogger("EDM4hepCalorimeterOutput", level)),
      m_config(config) {
  const std::set<std::string> names{config.inputHits, config.inputDeposits,
                                    config.inputSources, config.outputHits,
                                    config.outputLinks};
  if (names.size() != 5 || names.contains("")) {
    throw std::invalid_argument(
        "Calorimeter converter names must be distinct and nonempty");
  }
  m_hits.initialize(config.inputHits);
  m_deposits.initialize(config.inputDeposits);
  m_sources.initialize(config.inputSources);
  m_outputHits.initialize(config.outputHits);
  m_outputLinks.initialize(config.outputLinks);
}

ProcessCode EDM4hepCalorimeterOutputConverter::execute(
    const AlgorithmContext& context) const {
  const auto& hits = m_hits(context);
  const auto& deposits = m_deposits(context);
  const auto& sources = m_sources(context);
  if (sources.size() != deposits.size()) {
    throw std::invalid_argument("Deposit/source collection sizes disagree");
  }
  // Validate the sidecar before using its handles to construct persistent
  // links.
  std::set<std::pair<std::uint32_t, int>> usedSources;
  for (std::size_t index = 0; index < deposits.size(); ++index) {
    const auto& deposit = deposits[index];
    const auto& source = sources[index];
    if (!source.hit.isAvailable() || !source.contribution.isAvailable() ||
        source.hit.getObjectID().index < 0 ||
        source.contribution.getObjectID().index < 0 ||
        deposit.cellId != source.hit.getCellID() ||
        !std::isfinite(deposit.depositedEnergy) ||
        deposit.depositedEnergy < 0 || !std::isfinite(deposit.time) ||
        deposit.depositedEnergy !=
            source.contribution.getEnergy() * Acts::UnitConstants::GeV ||
        deposit.time !=
            source.contribution.getTime() * Acts::UnitConstants::ns ||
        !deposit.position.allFinite()) {
      throw std::invalid_argument("Deposit/source mapping is inconsistent");
    }
    const auto contributions = source.hit.getContributions();
    const auto objectId = source.contribution.getObjectID();
    if (std::ranges::find(contributions, source.contribution) ==
            contributions.end() ||
        !usedSources.emplace(objectId.collectionID, objectId.index).second) {
      throw std::invalid_argument(
          "Calorimeter source is repeated or belongs to another hit");
    }
  }
  edm4hep::CalorimeterHitCollection outputHits;
  edm4hep::CaloHitSimCaloHitLinkCollection outputLinks;
  std::set<std::uint64_t> usedCells;
  std::set<std::size_t> usedDeposits;
  for (const auto& hit : hits) {
    if (!std::isfinite(hit.energy) || hit.energy <= 0 ||
        !std::isfinite(hit.time) || !hit.position.allFinite() ||
        hit.sourceIndices.empty() || !usedCells.insert(hit.cellId).second) {
      throw std::invalid_argument("Invalid calibrated calorimeter hit");
    }
    double totalEnergy = 0;
    std::map<edm4hep::SimCalorimeterHit, double> contributions;
    for (const auto index : hit.sourceIndices) {
      if (index >= deposits.size() || !usedDeposits.insert(index).second ||
          deposits[index].cellId != hit.cellId ||
          deposits[index].depositedEnergy <= 0) {
        throw std::invalid_argument(
            "Invalid calibrated-hit deposit provenance");
      }
      contributions[sources[index].hit] += deposits[index].depositedEnergy;
      totalEnergy += deposits[index].depositedEnergy;
    }
    if (!std::isfinite(totalEnergy)) {
      throw std::overflow_error("Calorimeter link energy sum overflow");
    }
    const float energy = checkedFloat(hit.energy / Acts::UnitConstants::GeV);
    if (energy == 0) {
      throw std::overflow_error(
          "Calorimeter energy underflows EDM4hep float range");
    }
    auto output = outputHits.create();
    output.setCellID(hit.cellId);
    output.setEnergy(energy);
    output.setTime(checkedFloat(hit.time / Acts::UnitConstants::ns));
    output.setPosition(
        {checkedFloat(hit.position.x() / Acts::UnitConstants::mm),
         checkedFloat(hit.position.y() / Acts::UnitConstants::mm),
         checkedFloat(hit.position.z() / Acts::UnitConstants::mm)});
    for (const auto& [sourceHit, energySum] : contributions) {
      auto link = outputLinks.create();
      link.setFrom(output);
      link.setTo(sourceHit);
      link.setWeight(checkedFloat(energySum / totalEnergy));
    }
  }
  m_outputHits(context, std::move(outputHits));
  m_outputLinks(context, std::move(outputLinks));
  return ProcessCode::SUCCESS;
}

std::vector<std::string> EDM4hepCalorimeterOutputConverter::collections()
    const {
  return {m_config.outputHits, m_config.outputLinks};
}

}  // namespace ActsExamples
