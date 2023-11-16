// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// DEBUG: THIS REQUIRES THE BINS TO BE SET TO phi:41 z:11

template <typename SpacePoint>
std::set<std::size_t> Acts::ATLASBottomBinFinder<SpacePoint>::findBins(
    std::size_t phiBin, std::size_t zBin,
    const Acts::SpacePointGrid<SpacePoint>* binnedSP) {
  std::set<std::size_t> neighbourBins =
      binnedSP->neighborHoodIndices({phiBin, zBin}, 1);
  if (zBin == 6) {
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin - 1, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin + 1, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin, zBin - 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin - 1, zBin - 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin + 1, zBin - 1}));
    return neighbourBins;
  }
  if (zBin > 6) {
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin - 1, zBin + 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin + 1, zBin + 1}));
  } else {
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin, zBin - 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin - 1, zBin - 1}));
    neighbourBins.erase(binnedSP->getGlobalBinIndex({phiBin + 1, zBin - 1}));
  }
  if (zBin == 4) {
    neighbourBins.insert(binnedSP->getGlobalBinIndex({phiBin, zBin + 2}));
    neighbourBins.insert(binnedSP->getGlobalBinIndex({phiBin - 1, zBin + 2}));
    neighbourBins.insert(binnedSP->getGlobalBinIndex({phiBin + 1, zBin + 2}));
  }
  if (zBin == 8) {
    neighbourBins.insert(binnedSP->getGlobalBinIndex({phiBin, zBin - 2}));
    neighbourBins.insert(binnedSP->getGlobalBinIndex({phiBin - 1, zBin - 2}));
    neighbourBins.insert(binnedSP->getGlobalBinIndex({phiBin + 1, zBin - 2}));
  }
  return neighbourBins;
}
