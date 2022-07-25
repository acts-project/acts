// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include <memory>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {

using NeighborhoodVector = boost::container::small_vector<size_t, 10>;

/// Iterates over the elements of all bins given
/// by the indices parameter in the given SpacePointGrid.
/// Fullfills the forward iterator.
class NeighborhoodIterator {
 public:
  using sp_it_t = typename std::vector<Acts::SpacePoint*>::const_iterator;

  NeighborhoodIterator() = delete;

  NeighborhoodIterator(NeighborhoodVector indices,
                       const SpacePointGrid* spgrid) {
    m_grid = spgrid;
    m_indices = indices;
    m_curInd = 0;
    if (m_indices.size() > m_curInd) {
      m_curIt = std::begin(spgrid->at(m_indices[m_curInd]));
      m_binEnd = std::end(spgrid->at(m_indices[m_curInd]));
    }
  }

  NeighborhoodIterator(NeighborhoodVector indices, const SpacePointGrid* spgrid,
                       size_t curInd, sp_it_t curIt) {
    m_grid = spgrid;
    m_indices = indices;
    m_curInd = curInd;
    m_curIt = curIt;
    if (m_indices.size() > m_curInd) {
      m_binEnd = std::end(spgrid->at(m_indices[m_curInd]));
    }
  }
  static NeighborhoodIterator begin(NeighborhoodVector indices,
                                    const SpacePointGrid* spgrid) {
    auto nIt = NeighborhoodIterator(indices, spgrid);
    // advance until first non-empty bin or last bin
    if (nIt.m_curIt == nIt.m_binEnd) {
      ++nIt;
    }
    return nIt;
  }

  NeighborhoodIterator(const NeighborhoodIterator& other) {
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
      if (m_curIt != m_binEnd) {
        return;
      }
    }
    // increase bin index m_curInd until you find non-empty bin
    // or until m_curInd >= m_indices.size()-1
    while (m_curIt == m_binEnd && m_indices.size() - 1 > m_curInd) {
      m_curInd++;
      m_curIt = std::begin(m_grid->at(m_indices[m_curInd]));
      m_binEnd = std::end(m_grid->at(m_indices[m_curInd]));
    }
  }

  Acts::SpacePoint* operator*() { return *m_curIt; }

  bool operator!=(const NeighborhoodIterator& other) {
    return m_curIt != other.m_curIt || m_curInd != other.m_curInd;
  }

  // iterators within current bin
  sp_it_t m_curIt;
  sp_it_t m_binEnd;
  // number of bins
  NeighborhoodVector m_indices;
  // current bin
  size_t m_curInd;
  const Acts::SpacePointGrid* m_grid;
};

/// @c Neighborhood Used to access iterators to access a group of bins
/// returned by a BinFinder.
/// Fulfills the range_expression interface
class Neighborhood {
 public:
  Neighborhood() = delete;
  Neighborhood(NeighborhoodVector indices, const SpacePointGrid* spgrid) {
    m_indices = indices;
    m_spgrid = spgrid;
  }
  NeighborhoodIterator begin() {
    return NeighborhoodIterator::begin(m_indices, m_spgrid);
  }
  NeighborhoodIterator end() {
    return NeighborhoodIterator(m_indices, m_spgrid, m_indices.size() - 1,
                                std::end(m_spgrid->at(m_indices.back())));
  }

 private:
  NeighborhoodVector m_indices;
  const SpacePointGrid* m_spgrid;
};

/// @c BinnedSPGroupIterator Allows to iterate over all groups of bins
/// a provided BinFinder can generate for each bin of a provided SPGrid
class BinnedSPGroupIterator {
 public:
  BinnedSPGroupIterator& operator++() {
    if (zIndex < phiZbins[1]) {
      zIndex++;
    } else {
      zIndex = 1;
      phiIndex++;
    }

    size_t this_zIndex = zIndex;
    if (not customZorder.empty())
      this_zIndex = customZorder.at(this_zIndex - 1);

    // set current & neighbor bins only if bin indices valid
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      currentBin = NeighborhoodVector{
          grid->globalBinFromLocalBins({phiIndex, this_zIndex})};
      bottomBinIndices =
          m_bottomBinFinder->findBins(phiIndex, this_zIndex, grid);
      topBinIndices = m_topBinFinder->findBins(phiIndex, this_zIndex, grid);
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

  bool operator!=(const BinnedSPGroupIterator& otherState) {
    return !(this->operator==(otherState));
  }

  Neighborhood middle() { return Neighborhood(currentBin, grid); }

  Neighborhood bottom() { return Neighborhood(bottomBinIndices, grid); }

  Neighborhood top() { return Neighborhood(topBinIndices, grid); }

  BinnedSPGroupIterator(const SpacePointGrid* spgrid, BinFinder* botBinFinder,
                        BinFinder* tBinFinder,
                        boost::container::small_vector<size_t, 20> bins = {}) {
    grid = spgrid;
    m_bottomBinFinder = botBinFinder;
    m_topBinFinder = tBinFinder;
    phiZbins = grid->numLocalBins();
    phiIndex = 1;
    zIndex = 1;
    customZorder = bins;
    // if m_bins vector was not defined, use z bin 1 (zIndex) to start the
    // iterator, otherwise use the first value in m_bins vector
    size_t this_zIndex = bins.empty() ? zIndex : bins.front();
    outputIndex = grid->globalBinFromLocalBins({phiIndex, this_zIndex});
    currentBin = NeighborhoodVector{
        grid->globalBinFromLocalBins({phiIndex, this_zIndex})};
    bottomBinIndices = m_bottomBinFinder->findBins(phiIndex, this_zIndex, grid);
    topBinIndices = m_topBinFinder->findBins(phiIndex, this_zIndex, grid);
  }

  BinnedSPGroupIterator(const SpacePointGrid* spgrid, BinFinder* botBinFinder,
                        BinFinder* tBinFinder, size_t phiInd, size_t zInd,
                        boost::container::small_vector<size_t, 20> bins = {}) {
    m_bottomBinFinder = botBinFinder;
    m_topBinFinder = tBinFinder;
    grid = spgrid;
    phiIndex = phiInd;
    zIndex = zInd;
    phiZbins = grid->numLocalBins();
    customZorder = bins;
    // if m_bins vector was not defined, use the next z bin (zInd), otherwise
    // use the z bin value stored in m_bins vector for a custom order
    size_t this_zIndex =
        bins.empty()
            ? zIndex
            : (zIndex <= phiZbins[1] ? bins.at(zIndex - 1) : bins.back());
    outputIndex = grid->globalBinFromLocalBins({phiIndex, this_zIndex});
    currentBin =
        NeighborhoodVector(grid->globalBinFromLocalBins({phiInd, this_zIndex}));
    if (phiIndex <= phiZbins[0] && zIndex <= phiZbins[1]) {
      bottomBinIndices =
          m_bottomBinFinder->findBins(phiIndex, this_zIndex, grid);
      topBinIndices = m_topBinFinder->findBins(phiIndex, this_zIndex, grid);
    }
  }

 private:
  // middle spacepoint bin
  NeighborhoodVector currentBin;
  NeighborhoodVector bottomBinIndices;
  NeighborhoodVector topBinIndices;
  const SpacePointGrid* grid;
  size_t phiIndex = 1;
  size_t zIndex = 1;
  size_t outputIndex = 0;
  std::array<long unsigned int, 2ul> phiZbins;
  BinFinder* m_bottomBinFinder;
  BinFinder* m_topBinFinder;
  boost::container::small_vector<size_t, 20> customZorder;
  // 	bool start = true;
};

/// @c BinnedSPGroup Provides access to begin and end BinnedSPGroupIterator
/// for given BinFinders and SpacePointGrid.
/// Fulfills the range_expression interface.
class BinnedSPGroup {
 public:
  BinnedSPGroup() = delete;

  template <typename spacepoint_iterator_t>
  BinnedSPGroup(spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
                std::shared_ptr<Acts::BinFinder> botBinFinder,
                std::shared_ptr<Acts::BinFinder> tBinFinder,
                std::unique_ptr<SpacePointGrid> grid,
                const SeedfinderConfig& _config);

  size_t size() { return m_binnedSP->size(); }

  BinnedSPGroupIterator begin() {
    return BinnedSPGroupIterator(m_binnedSP.get(), m_bottomBinFinder.get(),
                                 m_topBinFinder.get(), m_bins);
  }

  BinnedSPGroupIterator end() {
    auto phiZbins = m_binnedSP->numLocalBins();
    return BinnedSPGroupIterator(m_binnedSP.get(), m_bottomBinFinder.get(),
                                 m_topBinFinder.get(), phiZbins[0],
                                 phiZbins[1] + 1, m_bins);
  }

 private:
  // grid holding pointers to all SpacePoints
  std::unique_ptr<Acts::SpacePointGrid> m_binnedSP;

  // BinFinder must return std::vector<Acts::Seeding::Bin> with content of
  // each bin sorted in r (ascending)
  std::shared_ptr<BinFinder> m_topBinFinder;
  std::shared_ptr<BinFinder> m_bottomBinFinder;

  boost::container::small_vector<size_t, 20> m_bins;
};

}  // namespace Acts
#include "Acts/Seeding/BinnedSPGroup.ipp"
