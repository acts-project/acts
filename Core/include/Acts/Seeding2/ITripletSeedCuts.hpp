// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"

namespace Acts::Experimental {

/// @c ITripletSeedCuts can be used to increase or decrease seed weights
/// based on the space points used in a seed. Seed weights are also
/// influenced by the SeedFilter default implementation. This tool is also used
/// to decide if a seed passes a seed weight cut. As the weight is stored in
/// seeds, there are two distinct methods.
class ITripletSeedCuts {
 public:
  virtual ~ITripletSeedCuts() = default;

  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  virtual float seedWeight(const ConstSpacePointProxy2& bottom,
                           const ConstSpacePointProxy2& middle,
                           const ConstSpacePointProxy2& top) const = 0;

  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  virtual bool singleSeedCut(float weight, const ConstSpacePointProxy2& bottom,
                             const ConstSpacePointProxy2& middle,
                             const ConstSpacePointProxy2& top) const = 0;

  /// @param seedCandidates contains collection of seed candidates created for one middle
  /// space point in a std::tuple format
  /// @return vector of seed candidates that pass the cut
  virtual void cutPerMiddleSp(
      std::span<TripletCandidate2> seedCandidates) const = 0;
};

}  // namespace Acts::Experimental
