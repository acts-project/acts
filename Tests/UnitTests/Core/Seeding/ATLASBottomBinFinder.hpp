// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/IBinFinder.hpp"

#include <set>

namespace Acts {

// DEBUG: THIS REQUIRES THE BINS TO BE SET TO phi:41 z:11

/// @class ATLASBinFinder
/// The ATLASBinFinder is an implementation of the ATLAS bincombination behavior
/// satisfying the IBinFinder interface. Assumes the grid has 11 bins filled by
/// the same logic as ATLAS bins.
template <typename SpacePoint>
class ATLASBottomBinFinder : public IBinFinder<SpacePoint> {
 public:
  /// destructor
  ~ATLASBottomBinFinder() = default;

  /// Return all bins that could contain space points that can be used with the
  /// space points in the bin with the provided indices to create seeds.
  /// @param phiBin phi index of bin with middle space points
  /// @param zBin z index of bin with middle space points
  /// @param binnedSP phi-z grid containing all bins
  std::set<std::size_t> findBins(std::size_t phiBin, std::size_t zBin,
                                 const SpacePointGrid<SpacePoint>* binnedSP);
};
}  // namespace Acts
#include "ATLASBottomBinFinder.ipp"
