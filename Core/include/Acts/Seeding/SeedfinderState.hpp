// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/IBinFinder.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include <memory>
#include <mutex>
#include <queue>
#include <vector>

namespace Acts {

template <typename SpacePoint>
class SeedfinderStateIterator
{
public:
  SeedfinderStateIterator& operator++()
  {
    if (zIndex < phiZbins[1]) {
      zIndex++;

    } else {
      zIndex = 1;
      phiIndex++;
    }
    // set current & neighbor bins only if bin indices valid
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      currentBin       = &(grid->at({phiIndex, zIndex}));
      bottomBinIndices = bottomBinFinder->findBins(phiIndex, zIndex, grid);
      topBinIndices    = topBinFinder->findBins(phiIndex, zIndex, grid);
      outputIndex++;
      return *this;
    }
    phiIndex = phiZbins[0];
    zIndex   = phiZbins[1] + 1;
    return *this;
  }

  bool
  operator==(const SeedfinderStateIterator& otherState)
  {
    return (zIndex == otherState.zIndex && phiIndex == otherState.phiIndex);
  }

  SeedfinderStateIterator(const SpacePointGrid<SpacePoint>* spgrid,
                       IBinFinder<SpacePoint>*           botBinFinder,
                       IBinFinder<SpacePoint>*           tBinFinder)
    : currentBin(&(spgrid->at({1, 1})))
  {
    grid             = spgrid;
    bottomBinFinder  = botBinFinder;
    topBinFinder     = tBinFinder;
    phiZbins         = grid->getNBins();
    phiIndex         = 1;
    zIndex           = 1;
    outputIndex      = 0;
    bottomBinIndices = bottomBinFinder->findBins(phiIndex, zIndex, grid);
    topBinIndices    = topBinFinder->findBins(phiIndex, zIndex, grid);
  }

  SeedfinderStateIterator(const SpacePointGrid<SpacePoint>* spgrid,
                       IBinFinder<SpacePoint>*           botBinFinder,
                       IBinFinder<SpacePoint>*           tBinFinder,
                       size_t                            phiInd,
                       size_t                            zInd)
    : currentBin(&(spgrid->at({phiInd, zInd})))
  {
    bottomBinFinder  = botBinFinder;
    topBinFinder     = tBinFinder;
    grid             = spgrid;
    phiIndex         = phiInd;
    zIndex           = zInd;
    phiZbins         = grid->getNBins();
    outputIndex      = (phiInd - 1) * phiZbins[1] + zInd - 1;
    bottomBinIndices = bottomBinFinder->findBins(phiIndex, zIndex, grid);
    topBinIndices    = topBinFinder->findBins(phiIndex, zIndex, grid);
  }

  // middle spacepoint bin
  const std::vector<std::unique_ptr<const InternalSpacePoint<SpacePoint>>>*
                                    currentBin;
  std::set<size_t>                  bottomBinIndices;
  std::set<size_t>                  topBinIndices;
  const SpacePointGrid<SpacePoint>* grid;
  size_t                            phiIndex    = 1;
  size_t                            zIndex      = 1;
  size_t                            outputIndex = 0;
  std::array<long unsigned int, 2ul> phiZbins;
  IBinFinder<SpacePoint>* bottomBinFinder;
  IBinFinder<SpacePoint>* topBinFinder;
};

template <typename SpacePoint>
struct SeedfinderState
{
  // grid with ownership of all InternalSpacePoint
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> binnedSP;

  // BinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<IBinFinder<SpacePoint>> bottomBinFinder;
  std::shared_ptr<IBinFinder<SpacePoint>> topBinFinder;

  // container with seeds created so far
  std::vector<std::vector<std::unique_ptr<Seed<SpacePoint>>>> outputVec;

  SeedfinderStateIterator<SpacePoint>
  begin()
  {
    return SeedfinderStateIterator<SpacePoint>(
        binnedSP.get(), bottomBinFinder.get(), topBinFinder.get());
  }

  SeedfinderStateIterator<SpacePoint>
  end()
  {
    auto phiZbins = binnedSP->getNBins();
    return SeedfinderStateIterator<SpacePoint>(binnedSP.get(),
                                            bottomBinFinder.get(),
                                            topBinFinder.get(),
                                            phiZbins[0],
                                            phiZbins[1] + 1);
  }
};
}
