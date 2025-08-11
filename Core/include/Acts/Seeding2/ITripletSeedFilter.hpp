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
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"

namespace Acts::Experimental {

class ITripletSeedFilter {
 public:
  virtual ~ITripletSeedFilter() = default;

  virtual void initialize(
      CandidatesForMiddleSp2& candidatesCollector) const = 0;

  virtual bool sufficientTopDoublets(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp& topDoublets) const = 0;

  /// Create seed candidates with fixed bottom and middle space points and
  /// all compatible top space points.
  ///
  /// @param spacePoints Container with all space points
  /// @param bottomLink Link to the bottom doublet space point
  /// @param spM Fixed middle space point
  /// @param topSpVec Vector containing all space points that may be compatible
  ///                 with both bottom and middle space point
  /// @param invHelixDiameterVec Vector containing 1/(2*r) values where r is the
  ///                            helix radius
  /// @param impactParametersVec Vector containing the impact parameters
  /// @param candidatesCollector Container for the seed candidates
  virtual void filterTripletTopCandidates(
      const SpacePointContainer2& spacePoints,
      const DoubletsForMiddleSp::Proxy& bottomLink,
      const ConstSpacePointProxy2& spM,
      std::span<const SpacePointIndex2> topSpVec,
      std::span<const float> invHelixDiameterVec,
      std::span<const float> impactParametersVec,
      CandidatesForMiddleSp2& candidatesCollector) const = 0;

  /// Create final seeds for all candidates with the same middle space point
  ///
  /// @param spacePoints Container with all space points
  /// @param candidates Collection of seed candidates
  /// @param numQualitySeeds Number of high quality seeds in seed confirmation
  /// @param outputCollection Output container for the seeds
  virtual void filterTripletsMiddleFixed(
      const SpacePointContainer2& spacePoints,
      std::span<TripletCandidate2> candidates, std::size_t numQualitySeeds,
      SeedContainer2& outputCollection) const = 0;
};

}  // namespace Acts::Experimental
