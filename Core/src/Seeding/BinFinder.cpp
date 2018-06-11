// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "Acts/Seeding/BinFinder.hpp"


std::set<size_t>
Acts::BinFinder::findBins(size_t i, size_t j, std::unique_ptr<Acts::SPGrid>& binnedSP){
  return binnedSP->neighborHoodIndices({i,j});  
}
