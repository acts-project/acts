// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SpacePointGrid.hpp"

namespace Acts {

/// @class IBinFinder
/// IBinFinder is the interface for both bottom bin finder as well as top bin
/// finder. It only contains one function that takes the indices of the bin
/// that contains all middle space points for the three space point seeds.
template <typename SpacePoint>
class IBinFinder
{
public:
  /// destructor
  virtual ~IBinFinder() = default;

  /// Returns all global bin indices (1D binning) that could contain bottom
  /// respectively top space points for the provided bin.
  /// @param phiBin the phi bin index of the bin containing the middle space
  /// points
  /// @param zBin the z bin index of the bin containing the middle space points
  /// @param binnedSP the grid containing all bins
  /// @return a set containing the global bin indices for all bins potentially
  /// containing bottom resp. top space points that can be combined with the
  /// middle space points from the provided bin
  virtual
  std::set<size_t>
  findBins(size_t                            phiBin,
           size_t                            zBin,
           const SpacePointGrid<SpacePoint>* binnedSP)
      = 0;
};
}
