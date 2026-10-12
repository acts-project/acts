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
#include "ActsCalorimeter/Geometry/CellNeighbours.hpp"

#include <map>
#include <span>

namespace ActsCalorimeter {

/// Deterministic seeded connected-component baseline, without shower splitting.
class CalorimeterClusterer {
 public:
  struct Config {
    /// Inclusive seed threshold on calibrated cell energy.
    double seedEnergyThreshold = 0;
    /// Inclusive neighbour threshold on calibrated cell energy. Zero-energy
    /// cells never join or bridge clusters, even with a zero threshold.
    double neighbourEnergyThreshold = 0;
  };

  /// Copy and normalize the fixed detector adjacency once. Repeated or reversed
  /// edges are harmless; self-edges are invalid. A cell absent from the graph
  /// is isolated. Geometry cells without an event hit do not bridge components.
  /// @throws std::invalid_argument for invalid thresholds or self-edges.
  explicit CalorimeterClusterer(const Config& config,
                                std::span<const CellNeighbour> neighbours = {});

  /// Connect positive cells at/above the neighbour threshold and retain each
  /// component containing at least one seed. Multiple seeds in a component
  /// produce one cluster. Each accepted cell belongs to exactly one cluster.
  /// Clusters are ordered by their smallest constituent cell ID. Energy sums
  /// and weighted means follow cell-ID order, independent of input/edge order.
  /// @throws std::invalid_argument for duplicate IDs or malformed cells.
  /// @throws std::overflow_error for nonfinite accumulated output.
  CalorimeterClusterContainer operator()(
      std::span<const CalorimeterHit> hits) const;

  const Config& config() const { return m_config; }

 private:
  Config m_config;
  std::map<std::uint64_t, std::vector<std::uint64_t>> m_neighbours;
};

}  // namespace ActsCalorimeter
