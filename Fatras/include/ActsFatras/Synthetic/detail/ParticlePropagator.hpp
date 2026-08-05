// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Walking the tracks of a synthetic event through a layout. See @ref
/// fatras_synthetic_events.

#include "Acts/Definitions/PdgParticle.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/detail/GeneratorScratch.hpp"

#include <random>
#include <vector>

namespace ActsFatras::Synthetic::detail {

/// Walks tracks through a layout: which surfaces each meets, what it measures
/// there, what the material does to it, and what it makes. Secondaries are
/// appended as they are made, so one call propagates the whole cascade.
class ParticlePropagator {
 public:
  /// @param layout the detector to propagate through, which has to outlive the
  ///        propagator
  /// @param config steering for the propagation
  /// @throws std::invalid_argument if `MaterialConfig::maxEnergyLossFraction`
  ///         is outside [0, 1), or if `EventConfig::particlePdg` names a
  ///         particle the core library has no mass for
  ParticlePropagator(const DetectorLayout& layout, const EventConfig& config);

  /// Propagate every track in `scratch.tracks`, appending the secondaries they
  /// make to it. `particles` holds what each track is recorded as, at the same
  /// index, and grows with it; only their hit counts are filled in here. The
  /// hits and the stubs of `scratch` are cleared first.
  ///
  /// @throws std::invalid_argument if `particles` holds fewer entries than
  ///         `scratch.tracks`
  void propagate(std::mt19937& rng, GeneratorScratch& scratch,
                 std::vector<GeneratedParticle>& particles) const;

 private:
  const DetectorLayout* m_layout{};
  EventConfig m_cfg;
  /// `EventConfig::particlePdg` without its sign, which is what the core
  /// library's energy loss is keyed on
  Acts::PdgParticle m_absPdg{};
  /// Mass of that species, which every track in the event is
  float m_mass{};
};

}  // namespace ActsFatras::Synthetic::detail
