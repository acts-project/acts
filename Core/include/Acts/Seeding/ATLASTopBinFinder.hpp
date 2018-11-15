// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include "Acts/Seeding/IBinFinder.hpp"

#include <set>

namespace Acts {

// DEBUG: THIS IS JUST A TEST TO SEE IF RESULT DIFFERENCES BETWEEN ACTS AND ATLAS ARE DUE TO HARDCODED BINFINDING


/// @class ATLASBinFinder
/// The ATLASBinFinder is an implementation of the ATLAS bincombination behavior
/// satisfying the IBinFinder interface. Assumes the grid has 11 bins filled by 
/// the same logic as ATLAS bins.
template <typename SpacePoint>
class ATLASTopBinFinder : public IBinFinder<SpacePoint>
{
public:
/// Virtual destructor
  ~ATLASTopBinFinder() = default;

/// Return all bins that could contain space points that can be used with the 
/// space points in the bin with the provided indices to create seeds.
/// @param phiBin phi index of bin with middle space points
/// @param zBin z index of bin with middle space points
/// @param binnedSP phi-z grid containing all bins
   virtual
   std::set<size_t>
   findBins(size_t phiBin, size_t zBin,const SpacePointGrid<SpacePoint>* binnedSP);

};
}

#include "Acts/Seeding/ATLASTopBinFinder.ipp"
