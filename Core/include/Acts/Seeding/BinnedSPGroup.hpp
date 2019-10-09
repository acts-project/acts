// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include <memory>
#include <vector>

namespace Acts {

template <typename external_spacepoint_t>
class NeighborhoodIterator {
 public:
  using sp_it_t = typename std::vector<std::unique_ptr<
      const InternalSpacePoint<external_spacepoint_t>>>::const_iterator;

  NeighborhoodIterator() = delete;

  NeighborhoodIterator(std::vector<size_t> indices,
                       const SpacePointGrid<external_spacepoint_t>* spgrid) {
    m_grid = spgrid;
    m_indices = indices;
    m_curInd = 0;
    if (m_indices.size() > m_curInd) {
      m_curIt = std::begin(spgrid->at(m_indices[m_curInd]));
      m_binEnd = std::end(spgrid->at(m_indices[m_curInd]));
    }
  }

  NeighborhoodIterator(std::vector<size_t> indices,
                       const SpacePointGrid<external_spacepoint_t>* spgrid,
                       size_t curInd, sp_it_t curIt) {
    m_grid = spgrid;
    m_indices = indices;
    m_curInd = curInd;
    m_curIt = curIt;
    if (m_indices.size() > m_curInd) {
      m_binEnd = std::end(spgrid->at(m_indices[m_curInd]));
    }
  }
  static NeighborhoodIterator<external_spacepoint_t> begin(std::vector<size_t> indices,                        
                                                           const SpacePointGrid<external_spacepoint_t>* spgrid) {
    auto nIt = NeighborhoodIterator<external_spacepoint_t>(indices, spgrid);
    // advance until first non-empty bin or last bin
    if(nIt.m_curIt == nIt.m_binEnd){
      ++nIt;
    }
    return nIt;
  }

  NeighborhoodIterator(
      const NeighborhoodIterator<external_spacepoint_t>& other) {
    m_grid = other.m_grid;
    m_indices = other.m_indices;
    m_curInd = other.m_curInd;
    m_curIt = other.m_curIt;
    m_binEnd = other.m_binEnd;
  }

  void operator++() {
    // if iterator of current Bin not yet at end, increase
    if (m_curIt != m_binEnd) {
      m_curIt++;
      // return only if end of current bin still not reached
      if (m_curIt != m_binEnd){
        return;
      }
    }
    // increase bin index m_curInd until you find non-empty bin
    // or until m_curInd >= m_indices.size()-1
    while(m_curIt == m_binEnd && m_indices.size()-1 > m_curInd){
      m_curInd++;
      m_curIt = std::begin(m_grid->at(m_indices[m_curInd]));
      m_binEnd = std::end(m_grid->at(m_indices[m_curInd]));
    } 
  }

  const InternalSpacePoint<external_spacepoint_t>* operator*() {
    return (*m_curIt).get();
  }

  bool operator==(const NeighborhoodIterator<external_spacepoint_t>& other) {
    return m_curIt == other.m_curIt && m_curInd == other.m_curInd;
  }

  // iterators within current bin
  sp_it_t m_curIt;
  sp_it_t m_binEnd;
  // number of bins
  std::vector<size_t> m_indices;
  // current bin
  size_t m_curInd;
  const Acts::SpacePointGrid<external_spacepoint_t>* m_grid;
};

template <typename external_spacepoint_t>
class BinnedSPGroupIterator {
 public:
  BinnedSPGroupIterator& operator++() {
    if (zIndex < phiZbins[1]) {
      zIndex++;

    } else {
      zIndex = 1;
      phiIndex++;
    }
    // set current & neighbor bins only if bin indices valid
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      currentBin =
          std::vector<size_t>{grid->globalBinFromLocalBins({phiIndex, zIndex})};
      bottomBinIndices = m_bottomBinFinder->findBins(phiIndex, zIndex, grid);
      topBinIndices = m_topBinFinder->findBins(phiIndex, zIndex, grid);
      outputIndex++;
      return *this;
    }
    phiIndex = phiZbins[0];
    zIndex = phiZbins[1] + 1;
    return *this;
  }

  bool operator==(const BinnedSPGroupIterator& otherState) {
    return (zIndex == otherState.zIndex && phiIndex == otherState.phiIndex);
  }

  NeighborhoodIterator<external_spacepoint_t> middleBegin() {
    return NeighborhoodIterator<external_spacepoint_t>::begin(currentBin, grid);
  }

  NeighborhoodIterator<external_spacepoint_t> middleEnd() {
    return NeighborhoodIterator<external_spacepoint_t>(
        currentBin, grid, currentBin.size()-1,
        std::end(grid->at(currentBin.back())));
  }

  NeighborhoodIterator<external_spacepoint_t> bottomBegin() {
    return NeighborhoodIterator<external_spacepoint_t>::begin(bottomBinIndices, grid);
  }

  NeighborhoodIterator<external_spacepoint_t> bottomEnd() {
    return NeighborhoodIterator<external_spacepoint_t>(
        bottomBinIndices, grid, bottomBinIndices.size()-1,
        std::end(grid->at(bottomBinIndices.back())));
  }

  NeighborhoodIterator<external_spacepoint_t> topBegin() {
    return NeighborhoodIterator<external_spacepoint_t>::begin(topBinIndices, grid);
  }

  NeighborhoodIterator<external_spacepoint_t> topEnd() {
    return NeighborhoodIterator<external_spacepoint_t>(
        topBinIndices, grid, topBinIndices.size()-1,
        std::end(grid->at(topBinIndices.back())));
  }

  BinnedSPGroupIterator(const SpacePointGrid<external_spacepoint_t>* spgrid,
                        BinFinder<external_spacepoint_t>* botBinFinder,
                        BinFinder<external_spacepoint_t>* tBinFinder)
      : currentBin({spgrid->globalBinFromLocalBins({1, 1})}) {
    grid = spgrid;
    m_bottomBinFinder = botBinFinder;
    m_topBinFinder = tBinFinder;
    phiZbins = grid->numLocalBins();
    phiIndex = 1;
    zIndex = 1;
    outputIndex = 0;
    bottomBinIndices = m_bottomBinFinder->findBins(phiIndex, zIndex, grid);
    topBinIndices = m_topBinFinder->findBins(phiIndex, zIndex, grid);
  }

  BinnedSPGroupIterator(const SpacePointGrid<external_spacepoint_t>* spgrid,
                        BinFinder<external_spacepoint_t>* botBinFinder,
                        BinFinder<external_spacepoint_t>* tBinFinder,
                        size_t phiInd, size_t zInd)
      : currentBin({spgrid->globalBinFromLocalBins({phiInd, zInd})}) {
    m_bottomBinFinder = botBinFinder;
    m_topBinFinder = tBinFinder;
    grid = spgrid;
    phiIndex = phiInd;
    zIndex = zInd;
    phiZbins = grid->numLocalBins();
    outputIndex = (phiInd - 1) * phiZbins[1] + zInd - 1;
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      bottomBinIndices = m_bottomBinFinder->findBins(phiIndex, zIndex, grid);
      topBinIndices = m_topBinFinder->findBins(phiIndex, zIndex, grid);
    }
  }

 private:
  // middle spacepoint bin
  std::vector<size_t> currentBin;
  std::vector<size_t> bottomBinIndices;
  std::vector<size_t> topBinIndices;
  const SpacePointGrid<external_spacepoint_t>* grid;
  size_t phiIndex = 1;
  size_t zIndex = 1;
  size_t outputIndex = 0;
  std::array<long unsigned int, 2ul> phiZbins;
  BinFinder<external_spacepoint_t>* m_bottomBinFinder;
  BinFinder<external_spacepoint_t>* m_topBinFinder;
};

template <typename external_spacepoint_t>
class BinnedSPGroup {
 public:
  BinnedSPGroup() = delete;

  template <typename spacepoint_iterator_t>
  BinnedSPGroup<external_spacepoint_t>(
      spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
      std::function<Acts::Vector2D(const external_spacepoint_t&, float, float,
                                   float)>
          covTool,
      std::shared_ptr<Acts::BinFinder<external_spacepoint_t>> botBinFinder,
      std::shared_ptr<Acts::BinFinder<external_spacepoint_t>> tBinFinder,
      std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
      SeedfinderConfig<external_spacepoint_t>& config);

  size_t size() { return m_binnedSP.size(); }

  BinnedSPGroupIterator<external_spacepoint_t> begin() {
    return BinnedSPGroupIterator<external_spacepoint_t>(
        m_binnedSP.get(), m_bottomBinFinder.get(), m_topBinFinder.get());
  }

  BinnedSPGroupIterator<external_spacepoint_t> end() {
    auto phiZbins = m_binnedSP->numLocalBins();
    return BinnedSPGroupIterator<external_spacepoint_t>(
        m_binnedSP.get(), m_bottomBinFinder.get(), m_topBinFinder.get(),
        phiZbins[0], phiZbins[1] + 1);
  }

 private:
  // grid with ownership of all InternalSpacePoint
  std::unique_ptr<Acts::SpacePointGrid<external_spacepoint_t>> m_binnedSP;

  // BinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<BinFinder<external_spacepoint_t>> m_topBinFinder;
  std::shared_ptr<BinFinder<external_spacepoint_t>> m_bottomBinFinder;
};

}  // namespace Acts
#include "Acts/Seeding/BinnedSPGroup.ipp"
