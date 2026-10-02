// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// The primaries of a synthetic pile-up event.

#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/detail/GeneratorScratch.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"

#include <random>
#include <vector>

namespace ActsFatras::Synthetic::detail {

/// Draws the primaries of one event: a transverse momentum, a rapidity and a
/// point in the luminous region each. It knows nothing of the detector.
class PrimaryGenerator {
 public:
  /// @param config steering for the generated primaries
  /// @throws std::invalid_argument if `EventConfig::particlePdg` names a
  ///         particle the core library has no mass for
  explicit PrimaryGenerator(const EventConfig& config);

  /// Append the primaries of one event to `tracks`, and what each of them is
  /// recorded as to `particles` at the same index.
  void generate(std::mt19937& rng, std::vector<Helix>& tracks,
                std::vector<GeneratedParticle>& particles) const;

 private:
  EventConfig m_cfg;
  /// Mass of `EventConfig::particlePdg`, the whole of the difference between a
  /// rapidity and a pseudorapidity
  float m_mass{};
};

}  // namespace ActsFatras::Synthetic::detail
