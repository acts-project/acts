// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SpacePointGrid.hpp"

#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @class BinFinder
/// The BinFinder is used by the ISPGroupSelector. It can be
/// used to find both bins that could be bottom bins as well as bins that could
/// be top bins, which are assumed to be the same bins. Does not take
/// interaction region into account to limit z-bins.
template <typename external_spacepoint_t>
class BinFinder {
 public:
  /// constructor
  BinFinder();

  BinFinder(const std::vector<std::pair<int, int> >&& zBinNeighbors,
            const int&& numPhiNeighbors);

  BinFinder(const std::vector<std::pair<int, int> >& zBinNeighbors,
            const int& numPhiNeighbors);

  /// destructor
  ~BinFinder() = default;

  /// Return all bins that could contain space points that can be used with the
  /// space points in the bin with the provided indices to create seeds.
  /// @param phiBin phi index of bin with middle space points
  /// @param zBin z index of bin with middle space points
  /// @param binnedSP phi-z grid containing all bins
  boost::container::small_vector<size_t, 10> findBins(
      size_t phiBin, size_t zBin,
      const SpacePointGrid<external_spacepoint_t>* binnedSP);

 private:
  const std::vector<std::pair<int, int> > m_zBinNeighbors;
  const int m_numPhiNeighbors;
};
}  // namespace Acts
#include "Acts/Seeding/BinFinder.ipp"
