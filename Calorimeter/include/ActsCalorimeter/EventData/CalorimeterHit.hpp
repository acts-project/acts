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

/// One simulated deposit. All quantities use ACTS native units.
struct SimCalorimeterHit {
  /// Full readout identifier; its bit encoding is owned by the detector.
  std::uint64_t cellId = 0;
  /// Cell centre in global coordinates, not the shower-step position.
  Acts::Vector3 position = Acts::Vector3::Zero();
  /// Deposited energy before response or shower-energy calibration.
  double depositedEnergy = 0;
  /// Deposit time in the input event's time convention; no TOF correction.
  double time = 0;
};

/// One calibrated cell after response and selection.
struct CalorimeterHit {
  std::uint64_t cellId = 0;
  Acts::Vector3 position = Acts::Vector3::Zero();
  /// Calibrated cell energy, distinct from simulated deposited energy.
  double energy = 0;
  /// Deposited-energy-weighted time of accepted deposits.
  double time = 0;
  /// Indices into the input SimCalorimeterHit collection for this event.
  /// These are transient provenance links, not persistent podio ObjectIDs.
  std::vector<std::size_t> sourceIndices;
};

using SimCalorimeterHitContainer = std::vector<SimCalorimeterHit>;
using CalorimeterHitContainer = std::vector<CalorimeterHit>;

}  // namespace ActsCalorimeter
