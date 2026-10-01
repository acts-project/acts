// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// A synthetic pile-up event on a simplified detector layout. See
/// @ref fatras_synthetic_events.

#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"
#include "ActsFatras/Synthetic/detail/GeneratorScratch.hpp"
#include "ActsFatras/Synthetic/detail/ParticlePropagator.hpp"
#include "ActsFatras/Synthetic/detail/PrimaryGenerator.hpp"

#include <cstdint>
#include <random>

namespace ActsFatras::Synthetic {

/// Generates a stream of synthetic events on one detector layout: successive
/// calls to `generate` continue the random sequence, `reset` restarts it.
///
/// Primaries are drawn, walked through the layout with everything they make,
/// and the hits they leave turned into space points.
class EventGenerator {
 public:
  /// @param layout the detector to generate hits on; not owned, it has to
  ///        outlive the generator
  /// @param config steering for the generated events
  /// @throws std::invalid_argument if `EventConfig::particlePdg` names a
  ///         particle the core library has no mass for, or if
  ///         `MaterialConfig::maxEnergyLossFraction` is outside [0, 1)
  EventGenerator(const DetectorLayout& layout, const EventConfig& config);

  /// Generate the next event, reusing the storage `event` already holds.
  /// @param event receives the generated event
  void generate(Event& event);

  /// Generate the next event.
  /// @return the generated event
  Event generate();

  /// Restart the random sequence, so that the next event is the first one
  /// again.
  ///
  /// @param seed the seed to restart from
  void reset(std::uint32_t seed);

 private:
  const DetectorLayout* m_layout{};
  EventConfig m_cfg;
  detail::PrimaryGenerator m_primaries;
  detail::ParticlePropagator m_propagator;
  std::mt19937 m_rng;
  detail::GeneratorScratch m_scratch;
};

/// Generate one event from a freshly seeded generator.
/// @param layout the detector to generate hits on
/// @param config steering for the generated event
/// @return the generated event
Event generateEvent(const DetectorLayout& layout, const EventConfig& config);

}  // namespace ActsFatras::Synthetic
