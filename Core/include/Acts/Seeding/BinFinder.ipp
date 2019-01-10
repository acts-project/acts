// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
template <typename SpacePoint>
std::set<size_t>
Acts::BinFinder<SpacePoint>::findBins(
    size_t                                  phiBin,
    size_t                                  zBin,
    const Acts::SpacePointGrid<SpacePoint>* binnedSP)
{
  return binnedSP->neighborHoodIndices({phiBin, zBin});
}
