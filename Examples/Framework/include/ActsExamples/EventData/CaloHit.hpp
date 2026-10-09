// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstdint>
#include <vector>

namespace ActsExamples {

/// Energy contribution to a calorimeter cell from a single source particle.
/// Already aggregated across all step-level depositions, so @c energy is the
/// total energy and @c time is the energy-weighted average.
struct CaloHitContribution {
  /// Row index of the source particle in the simulated
  /// @c SimParticleContainer; @c std::numeric_limits<std::uint64_t>::max() if
  /// the source particle was filtered out upstream.
  std::uint64_t particleRow;
  /// Sum of step-level energies in the same units as the source EDM4hep
  /// (typically GeV).
  float energy;
  /// Energy-weighted average step time in ns.
  float time;
};

/// A single calorimeter cell with its per-particle contribution list. Cells
/// are produced by an upstream input converter that has already applied any
/// time-of-flight correction and timing window cuts; downstream writers may
/// still drop cells based on energy thresholds.
struct CaloHit {
  /// Sensor cell identifier (the EDM4hep @c cellID).
  std::uint64_t cellId{};
  /// Subdetector enum chosen by the input converter (e.g. ECal/HCal × side).
  /// 255 is reserved for "unknown".
  std::uint8_t detector{};
  /// Cell position in mm.
  Acts::Vector3 position{};
  /// Sum of all surviving contributions' energies.
  float totalEnergy{};
  /// Per-particle contributions, in stable insertion order.
  std::vector<CaloHitContribution> contributions;
};

using CaloHitContainer = std::vector<CaloHit>;

}  // namespace ActsExamples
