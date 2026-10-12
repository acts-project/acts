// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsCalorimeter/Digitization/CalorimeterResponse.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace ActsCalorimeter {

CalorimeterResponse::CalorimeterResponse(const Config& config)
    : m_config(config) {
  if (!std::isfinite(config.energyScale) || config.energyScale <= 0 ||
      !std::isfinite(config.energyThreshold) || config.energyThreshold < 0 ||
      std::isnan(config.timeMin) || std::isnan(config.timeMax) ||
      config.timeMin > config.timeMax ||
      config.timeMin == std::numeric_limits<double>::infinity() ||
      config.timeMax == -std::numeric_limits<double>::infinity()) {
    throw std::invalid_argument("Invalid calorimeter response configuration");
  }
}

CalorimeterHitContainer CalorimeterResponse::operator()(
    std::span<const SimCalorimeterHit> deposits) const {
  std::map<std::uint64_t, CalorimeterHit> cells;
  for (std::size_t index = 0; index < deposits.size(); ++index) {
    const auto& deposit = deposits[index];
    if (!std::isfinite(deposit.depositedEnergy) ||
        deposit.depositedEnergy < 0 || !std::isfinite(deposit.time) ||
        !deposit.position.allFinite()) {
      throw std::invalid_argument("Invalid simulated calorimeter deposit");
    }
    if (deposit.depositedEnergy == 0 || deposit.time < m_config.timeMin ||
        deposit.time > m_config.timeMax) {
      continue;
    }

    auto [entry, inserted] = cells.try_emplace(deposit.cellId);
    auto& cell = entry->second;
    if (inserted) {
      cell.cellId = deposit.cellId;
      cell.position = deposit.position;
      cell.time = deposit.time;
    } else if (!(cell.position.array() == deposit.position.array()).all()) {
      throw std::invalid_argument("Inconsistent centre for calorimeter cell");
    }
    // During accumulation this holds deposited energy; calibration follows.
    const double total = cell.energy + deposit.depositedEnergy;
    if (!std::isfinite(total)) {
      throw std::overflow_error("Calorimeter deposited-energy sum overflow");
    }
    const double fraction = deposit.depositedEnergy / total;
    cell.time = (1 - fraction) * cell.time + fraction * deposit.time;
    cell.energy = total;
    cell.sourceIndices.push_back(index);
  }

  CalorimeterHitContainer output;
  output.reserve(cells.size());
  for (auto& [cellId, cell] : cells) {
    cell.energy *= m_config.energyScale;
    if (!std::isfinite(cell.energy) || !std::isfinite(cell.time)) {
      throw std::overflow_error("Calorimeter calibrated output overflow");
    }
    if (cell.energy > 0 && cell.energy >= m_config.energyThreshold) {
      output.push_back(std::move(cell));
    }
  }
  return output;
}

}  // namespace ActsCalorimeter
