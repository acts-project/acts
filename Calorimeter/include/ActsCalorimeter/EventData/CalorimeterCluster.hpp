// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace ActsCalorimeter {

/// A connected component of calibrated calorimeter cells in ACTS native units.
struct CalorimeterCluster {
  /// Sum of constituent calibrated energies; no further calibration is applied.
  double energy = 0;
  /// Calibrated-energy-weighted global cell centre.
  Acts::Vector3 position = Acts::Vector3::Zero();
  /// Calibrated-energy-weighted cell time, without TOF correction.
  double time = 0;
  /// Highest-energy cell; ties select the smaller cell ID.
  std::uint64_t seedCellId = 0;
  /// Indices into this event's input CalorimeterHit collection, sorted by cell
  /// ID. Follow each hit's sourceIndices to reach simulated deposits. These
  /// transient indices are not persistent podio ObjectIDs.
  std::vector<std::size_t> hitIndices;
};

using CalorimeterClusterContainer = std::vector<CalorimeterCluster>;

}  // namespace ActsCalorimeter
