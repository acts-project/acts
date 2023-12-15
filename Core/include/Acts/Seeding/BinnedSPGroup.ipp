// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Binned SP Group Iterator

#include <boost/container/flat_set.hpp>

template <typename external_spacepoint_t>
Acts::BinnedSPGroupIterator<external_spacepoint_t>::BinnedSPGroupIterator(
    Acts::BinnedSPGroup<external_spacepoint_t>& group,
    std::array<std::size_t, 2> index,
    std::array<std::vector<std::size_t>, 2> navigation)
    : m_group(group), m_gridItr(*group.m_grid.get(), index, navigation) {
  std::array<std::size_t, 2ul> endline{};
  endline[0ul] = navigation[0ul].size();
  endline[1ul] = navigation[1ul].size();
  m_gridItrEnd =
      typename Acts::SpacePointGrid<external_spacepoint_t>::local_iterator_t(
          *group.m_grid.get(), endline, std::move(navigation));
  findNotEmptyBin();
}

template <typename external_spacepoint_t>
inline Acts::BinnedSPGroupIterator<external_spacepoint_t>&
Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator++() {
  ++m_gridItr;
  findNotEmptyBin();
  return *this;
}

template <typename external_spacepoint_t>
inline bool Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator==(
    const Acts::BinnedSPGroupIterator<external_spacepoint_t>& other) const {
  return m_group.ptr == other.m_group.ptr && m_gridItr == other.m_gridItr;
}

template <typename external_spacepoint_t>
inline bool Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator!=(
    const Acts::BinnedSPGroupIterator<external_spacepoint_t>& other) const {
  return !(*this == other);
}

template <typename external_spacepoint_t>
std::tuple<boost::container::small_vector<std::size_t, 9>, std::size_t,
           boost::container::small_vector<std::size_t, 9>>
Acts::BinnedSPGroupIterator<external_spacepoint_t>::operator*() const {
  // Global Index
  std::array<std::size_t, 2> localPosition = m_gridItr.localPosition();
  std::size_t global_index =
      m_group->m_grid->globalBinFromLocalBins(localPosition);

  auto bottoms =
      m_group->m_bottomBinFinder->findBins(localPosition,
                                           *m_group->m_grid.get());
  auto tops =
      m_group->m_topBinFinder->findBins(localPosition,
                                        *m_group->m_grid.get());

  // GCC12+ in Release throws an overread warning here due to the move.
  // This is from inside boost code, so best we can do is to suppress it.
#if defined(__GNUC__) && __GNUC__ >= 12 && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  return std::make_tuple(std::move(bottoms), global_index, std::move(tops));
#if defined(__GNUC__) && __GNUC__ >= 12 && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
}

template <typename external_spacepoint_t>
inline void
Acts::BinnedSPGroupIterator<external_spacepoint_t>::findNotEmptyBin() {
  if (m_gridItr == m_gridItrEnd) {
    return;
  }
  // Iterate on the grid till we find a not-empty bin
  // We start from the current bin configuration and move forward
  std::size_t dimCollection = (*m_gridItr).size();
  while (dimCollection == 0ul && ++m_gridItr != m_gridItrEnd) {
    dimCollection = (*m_gridItr).size();
  }
}

// Binned SP Group
template <typename external_spacepoint_t>
template <typename spacepoint_iterator_t, typename callable_t>
Acts::BinnedSPGroup<external_spacepoint_t>::BinnedSPGroup(
    spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
    callable_t&& toGlobal,
    std::shared_ptr<const Acts::GridBinFinder<2ul>> botBinFinder,
    std::shared_ptr<const Acts::GridBinFinder<2ul>> tBinFinder,
    std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
    Acts::Extent& rRangeSPExtent,
    const SeedFinderConfig<external_spacepoint_t>& config,
    const SeedFinderOptions& options) {
  if (!config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderConfig not in ACTS internal units in BinnedSPGroup");
  }
  if (!options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptions not in ACTS internal units in BinnedSPGroup");
  }
  static_assert(
      std::is_same<
          typename std::iterator_traits<spacepoint_iterator_t>::value_type,
          const external_spacepoint_t*>::value,
      "Iterator does not contain type this class was templated with");

  // get region of interest (or full detector if configured accordingly)
  float phiMin = config.phiMin;
  float phiMax = config.phiMax;
  float zMin = config.zMin;
  float zMax = config.zMax;

  // sort by radius
  // add magnitude of beamPos to rMax to avoid excluding measurements
  // create number of bins equal to number of millimeters rMax
  // (worst case minR: configured minR + 1mm)
  // binSizeR allows to increase or reduce numRBins if needed
  std::size_t numRBins = static_cast<std::size_t>(
      (config.rMax + options.beamPos.norm()) / config.binSizeR);

  // keep track of changed bins while sorting
  boost::container::flat_set<std::size_t> rBinsIndex;

  std::size_t counter = 0;
  for (spacepoint_iterator_t it = spBegin; it != spEnd; it++, ++counter) {
    if (*it == nullptr) {
      continue;
    }
    const external_spacepoint_t& sp = **it;
    const auto& [spPosition, variance] =
        toGlobal(sp, config.zAlign, config.rAlign, config.sigmaError);

    float spX = spPosition[0];
    float spY = spPosition[1];
    float spZ = spPosition[2];

    // store x,y,z values in extent
    rRangeSPExtent.extend({spX, spY, spZ});

    // remove SPs outside z and phi region
    if (spZ > zMax || spZ < zMin) {
      continue;
    }
    float spPhi = std::atan2(spY, spX);
    if (spPhi > phiMax || spPhi < phiMin) {
      continue;
    }

    auto isp = std::make_unique<InternalSpacePoint<external_spacepoint_t>>(
        counter, sp, spPosition, options.beamPos, variance);
    // calculate r-Bin index and protect against overflow (underflow not
    // possible)
    std::size_t rIndex =
        static_cast<std::size_t>(isp->radius() / config.binSizeR);
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins) {
      continue;
    }

    // fill rbins into grid
    Acts::Vector2 spLocation(isp->phi(), isp->z());
    std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>&
        rbin = grid->atPosition(spLocation);
    rbin.push_back(std::move(isp));

    // keep track of the bins we modify so that we can later sort the SPs in
    // those bins only
    if (rbin.size() > 1) {
      rBinsIndex.insert(grid->globalBinFromPosition(spLocation));
    }
  }

  // sort SPs in R for each filled (z, phi) bin
  for (auto& binIndex : rBinsIndex) {
    std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>&
        rbin = grid->atPosition(binIndex);
    std::sort(
        rbin.begin(), rbin.end(),
        [](std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>& a,
           std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>& b) {
          return a->radius() < b->radius();
        });
  }

  m_grid = std::move(grid);
  m_bottomBinFinder = botBinFinder;
  m_topBinFinder = tBinFinder;

  m_skipZMiddleBin = config.skipZMiddleBinSearch;

  // phi axis
  m_bins[INDEX::PHI].resize(m_grid->numLocalBins()[0]);
  std::iota(m_bins[INDEX::PHI].begin(), m_bins[INDEX::PHI].end(), 1ul);

  // z axis
  if (config.zBinsCustomLooping.empty()) {
    std::size_t nZbins = m_grid->numLocalBins()[1] - m_skipZMiddleBin;
    m_bins[INDEX::Z] = std::vector<std::size_t>(nZbins);
    std::iota(m_bins[INDEX::Z].begin(), m_bins[INDEX::Z].end(),
              1ul + m_skipZMiddleBin);
  } else {
    m_bins[INDEX::Z] = std::vector<std::size_t>(
        config.zBinsCustomLooping.begin() + m_skipZMiddleBin,
        config.zBinsCustomLooping.end());
  }
}

template <typename external_spacepoint_t>
inline std::size_t Acts::BinnedSPGroup<external_spacepoint_t>::size() const {
  return m_grid->size();
}

template <typename external_spacepoint_t>
inline Acts::BinnedSPGroupIterator<external_spacepoint_t>
Acts::BinnedSPGroup<external_spacepoint_t>::begin() {
  return {*this, {0ul, 0ul}, m_bins};
}

template <typename external_spacepoint_t>
inline Acts::BinnedSPGroupIterator<external_spacepoint_t>
Acts::BinnedSPGroup<external_spacepoint_t>::end() {
  std::array<std::size_t, 2ul> endline{};
  endline[0ul] = m_bins[0ul].size();
  endline[1ul] = m_bins[1ul].size();
  return {*this, endline, m_bins};
}
