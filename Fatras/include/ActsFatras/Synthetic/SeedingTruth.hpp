// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file
/// Truth bookkeeping for a synthetic event, cheap enough to run alongside a
/// benchmark.

#include "Acts/EventData/SeedContainer.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"

#include <cstddef>

namespace ActsFatras::Synthetic {

/// Truth content of an event, so that a caller can report what it is running
/// on.
struct EventSummary {
  /// Number of space points
  std::size_t spacePoints{};
  /// Number of generated primaries
  std::size_t primaries{};
  /// Number of generated secondaries
  std::size_t secondaries{};
  /// Primaries above the momentum threshold that leave enough space points to
  /// be seedable at all
  std::size_t seedablePrimaries{};
  /// Space points left by primaries
  std::size_t primaryHits{};
  /// Space points left by secondaries
  std::size_t secondaryHits{};
};

/// Summarise the truth content of an event.
/// @param event the generated event
/// @param ptThreshold the momentum threshold a primary must pass
/// @return the summary
EventSummary summarize(const Event& event, float ptThreshold);

/// What the seeds of an event look like against the generator truth.
struct SeedingSummary {
  /// Number of seeds
  std::size_t seeds{};
  /// Seeds sharing enough space points with one primary above the momentum
  /// threshold
  std::size_t trueSeeds{};
  /// Primaries above the momentum threshold that have at least one true seed
  std::size_t matchedPrimaries{};
};

/// Match seeds against the generator truth.
/// @param event the generated event
/// @param seeds the seeds found on it, indexing `event.spacePoints`
/// @param ptThreshold the momentum threshold a primary must pass
/// @param minTrueSpacePoints space points a seed has to share with one particle
///        to be its seed; three constrain a helix
/// @return the summary
SeedingSummary evaluateSeeds(const Event& event,
                             const Acts::SeedContainer& seeds,
                             float ptThreshold,
                             std::size_t minTrueSpacePoints = 3);

}  // namespace ActsFatras::Synthetic
