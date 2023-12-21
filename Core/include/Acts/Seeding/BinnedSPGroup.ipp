// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Binned SP Group Iterator

#include <boost/container/flat_set.hpp>

namespace Acts {

template <typename grid_t>
BinnedGroup<grid_t>::BinnedGroup(std::unique_ptr<grid_t> grid,
				 std::shared_ptr<const Acts::GridBinFinder<Acts::BinnedGroup<grid_t>::DIM>> bottomFinder,
				 std::shared_ptr<const Acts::GridBinFinder<Acts::BinnedGroup<grid_t>::DIM>> topFinder,
				 std::array<std::vector<std::size_t>, Acts::BinnedGroup<grid_t>::DIM> navigation)
  : m_grid(std::move(grid)),
    m_bottomBinFinder(std::move(bottomFinder)),
    m_topBinFinder(std::move(topFinder)),
    m_bins(std::move(navigation))
{
  assert(m_grid != nullptr);
  assert(m_bottomBinFinder != nullptr);
  assert(m_topBinFinder != nullptr);

  /// If navigation is now defined for all axes, then we default that to a std::iota from 1ul
  std::array<std::size_t, DIM> numLocBins = m_grid->numLocalBins();
  for (std::size_t i(0ul); i<DIM; ++i) {
    if (!m_bins[i].empty()) {
      continue;
    }
    m_bins[i].resize(numLocBins[i]);
    std::iota(m_bins[i].begin(), m_bins[i].end(), 1ul);
  }
}
  
template <typename grid_t>
const grid_t& BinnedGroup<grid_t>::grid() const { return *m_grid.get(); }

template <typename grid_t>
grid_t& BinnedGroup<grid_t>::grid() { return *m_grid.get(); }

template <typename grid_t>
Acts::BinnedGroupIterator<grid_t> BinnedGroup<grid_t>::begin() const {
  return Acts::BinnedGroupIterator<grid_t>(*this, std::array<std::size_t, Acts::BinnedGroup<grid_t>::DIM>(), m_bins);
}

template <typename grid_t>
Acts::BinnedGroupIterator<grid_t> BinnedGroup<grid_t>::end() const {
  std::array<std::size_t, Acts::BinnedGroup<grid_t>::DIM> endline{};
  for (std::size_t i(0ul); i<Acts::BinnedGroup<grid_t>::DIM; ++i) {
    endline[i] = m_bins[i].size();
  }
  return Acts::BinnedGroupIterator<grid_t>(*this, endline, m_bins);
}

template <typename grid_t>
template <typename external_spacepoint_t,
	  typename external_spacepoint_iterator_t,
	  typename callable_t>
void BinnedGroup<grid_t>::fill(const Acts::SeedFinderConfig<external_spacepoint_t>& config,
			       const Acts::SeedFinderOptions& options,
			       external_spacepoint_iterator_t spBegin, external_spacepoint_iterator_t spEnd,
			       callable_t&& toGlobal,
			       Acts::Extent& rRangeSPExtent) {

  using iterated_value_t = typename std::iterator_traits<external_spacepoint_iterator_t>::value_type;
  using iterated_t = typename std::remove_const<typename std::remove_pointer<iterated_value_t>::type>::type;
  static_assert(std::is_same<iterated_t, external_spacepoint_t>::value,
		"Iterator does not contain type this class was templated with");
  
  if (!config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderConfig not in ACTS internal units in BinnedSPGroup");
  }
  if (!options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptions not in ACTS internal units in BinnedSPGroup");
  }

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

  std::size_t counter = 0ul;
  for (external_spacepoint_iterator_t it = spBegin; it != spEnd; it++, ++counter) {
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
        rbin = m_grid->atPosition(spLocation);
    rbin.push_back(std::move(isp));

    // keep track of the bins we modify so that we can later sort the SPs in
    // those bins only
    if (rbin.size() > 1) {
      rBinsIndex.insert(m_grid->globalBinFromPosition(spLocation));
    }
  }
  
  /// sort SPs in R for each filled bin
  for (auto& binIndex : rBinsIndex) {
    auto& rbin = m_grid->atPosition(binIndex);
    std::sort(
        rbin.begin(), rbin.end(),
        [](const auto& a, const auto& b) {
          return a->radius() < b->radius();
        });
  }  
}

} // namespace Acts
