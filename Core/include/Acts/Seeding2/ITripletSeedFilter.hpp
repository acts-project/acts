// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"

namespace Acts {

/// Interface for triplet seed filtering.
///
/// The filter is expected to be stateful and maintain internal information
/// across calls.
class ITripletSeedFilter {
 public:
  virtual ~ITripletSeedFilter() = default;

  /// @brief Check if there are sufficient top doublets for triplet formation
  /// @param spacePoints Container of space points
  /// @param spM Middle space point proxy
  /// @param topDoublets Collection of top doublets for the middle space point
  /// @return True if sufficient doublets are available for triplet seeds
  virtual bool sufficientTopDoublets(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp& topDoublets) const = 0;

  /// Create seed candidates with fixed bottom and middle space points and
  /// all compatible top space points.
  ///
  /// @param spacePoints Container with all space points
  /// @param spM Fixed middle space point
  /// @param bottomLink Link to the bottom doublet space point
  /// @param tripletTopCandidates Collection of triplet top candidates
  virtual void filterTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomLink,
      const TripletTopCandidates& tripletTopCandidates) const = 0;

  /// Create final seeds for all candidates with the same middle space point
  ///
  /// @param spacePoints Container with all space points
  /// @param outputCollection Output container for the seeds
  virtual void filterTripletsMiddleFixed(
      const SpacePointContainer2& spacePoints,
      SeedContainer2& outputCollection) const = 0;
};

}  // namespace Acts
