// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CandidatesForMiddleSp.hpp"
#include "Acts/Seeding/InternalSpacePointContainer.hpp"

namespace Acts {

/// @c IExperimentCuts can be used to increase or decrease seed weights
/// based on the space points used in a seed. Seed weights are also
/// influenced by the SeedFilter default implementation. This tool is also used
/// to decide if a seed passes a seed weight cut. As the weight is stored in
/// seeds, there are two distinct methods.
class IExperimentCuts {
 public:
  virtual ~IExperimentCuts() = default;

  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  virtual float seedWeight(const ConstInternalSpacePointProxy& bottom,
                           const ConstInternalSpacePointProxy& middle,
                           const ConstInternalSpacePointProxy& top) const = 0;
  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  virtual bool singleSeedCut(float weight,
                             const ConstInternalSpacePointProxy& bottom,
                             const ConstInternalSpacePointProxy& middle,
                             const ConstInternalSpacePointProxy& top) const = 0;

  /// @param seedCandidates contains collection of seed candidates created for one middle
  /// space point in a std::tuple format
  /// @param spacePoints container for the space points
  /// @return vector of seed candidates that pass the cut
  virtual std::vector<TripletCandidate> cutPerMiddleSP(
      std::vector<TripletCandidate> seedCandidates,
      const InternalSpacePointContainer& spacePoints) const = 0;
};

}  // namespace Acts
