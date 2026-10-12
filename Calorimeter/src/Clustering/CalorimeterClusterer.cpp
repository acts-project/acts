// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsCalorimeter/Clustering/CalorimeterClusterer.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace ActsCalorimeter {

CalorimeterClusterer::CalorimeterClusterer(
    const Config& config, std::span<const CellNeighbour> neighbours)
    : m_config(config) {
  if (!std::isfinite(config.seedEnergyThreshold) ||
      !std::isfinite(config.neighbourEnergyThreshold) ||
      config.neighbourEnergyThreshold < 0 ||
      config.seedEnergyThreshold < config.neighbourEnergyThreshold) {
    throw std::invalid_argument("Invalid calorimeter clustering thresholds");
  }
  for (const auto& edge : neighbours) {
    if (edge.first == edge.second) {
      throw std::invalid_argument("Calorimeter cell cannot neighbour itself");
    }
    m_neighbours[edge.first].push_back(edge.second);
    m_neighbours[edge.second].push_back(edge.first);
  }
  for (auto& [cellId, adjacent] : m_neighbours) {
    std::ranges::sort(adjacent);
    adjacent.erase(std::unique(adjacent.begin(), adjacent.end()),
                   adjacent.end());
  }
}

CalorimeterClusterContainer CalorimeterClusterer::operator()(
    std::span<const CalorimeterHit> hits) const {
  std::map<std::uint64_t, std::size_t> cells;
  for (std::size_t index = 0; index < hits.size(); ++index) {
    const auto& hit = hits[index];
    if (!std::isfinite(hit.energy) || hit.energy < 0 ||
        !std::isfinite(hit.time) || !hit.position.allFinite()) {
      throw std::invalid_argument("Invalid calibrated calorimeter cell");
    }
    if (!cells.emplace(hit.cellId, index).second) {
      throw std::invalid_argument("Duplicate calibrated calorimeter cell ID");
    }
  }
  auto eligible = [&](std::size_t index) {
    return hits[index].energy > 0 &&
           hits[index].energy >= m_config.neighbourEnergyThreshold;
  };

  CalorimeterClusterContainer output;
  std::vector<bool> visited(hits.size(), false);
  for (const auto& [cellId, index] : cells) {
    if (visited[index] || !eligible(index)) {
      continue;
    }
    // The growing component doubles as an iterative traversal queue.
    std::vector<std::size_t> members{index};
    visited[index] = true;
    for (std::size_t next = 0; next < members.size(); ++next) {
      const auto neighbours = m_neighbours.find(hits[members[next]].cellId);
      if (neighbours == m_neighbours.end()) {
        continue;
      }
      for (const auto neighbourId : neighbours->second) {
        const auto found = cells.find(neighbourId);
        if (found != cells.end() && !visited[found->second] &&
            eligible(found->second)) {
          visited[found->second] = true;
          members.push_back(found->second);
        }
      }
    }
    if (std::ranges::none_of(members, [&](std::size_t member) {
          return hits[member].energy >= m_config.seedEnergyThreshold;
        })) {
      continue;
    }
    std::ranges::sort(members, [&](std::size_t left, std::size_t right) {
      return hits[left].cellId < hits[right].cellId;
    });

    CalorimeterCluster cluster;
    double highestEnergy = -1;
    for (const auto member : members) {
      const auto& hit = hits[member];
      const double total = cluster.energy + hit.energy;
      if (!std::isfinite(total)) {
        throw std::overflow_error("Calorimeter cluster energy sum overflow");
      }
      const double fraction = hit.energy / total;
      cluster.position =
          (1 - fraction) * cluster.position + fraction * hit.position;
      cluster.time = (1 - fraction) * cluster.time + fraction * hit.time;
      cluster.energy = total;
      if (hit.energy > highestEnergy) {
        highestEnergy = hit.energy;
        cluster.seedCellId = hit.cellId;
      }
    }
    if (!cluster.position.allFinite() || !std::isfinite(cluster.time)) {
      throw std::overflow_error("Calorimeter cluster weighted-mean overflow");
    }
    cluster.hitIndices = std::move(members);
    output.push_back(std::move(cluster));
  }
  return output;
}

}  // namespace ActsCalorimeter
