// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <chrono>

template <typename external_spacepoint_t>
template <typename container_t, typename callable_t>
Acts::BinnedSPGroup<external_spacepoint_t>::BinnedSPGroup(
    const container_t& spacePoints, callable_t&& toGlobal,
    std::shared_ptr<const Acts::BinFinder<external_spacepoint_t>> botBinFinder,
    std::shared_ptr<const Acts::BinFinder<external_spacepoint_t>> tBinFinder,
    std::unique_ptr<SpacePointGrid<external_spacepoint_t>> grid,
    Acts::Extent& rRangeSPExtent,
    const SeedFinderConfig<external_spacepoint_t>& config,
    const SeedFinderOptions& options) {
  if (not config.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderConfig not in ACTS internal units in BinnedSPGroup");
  }
  if (not options.isInInternalUnits) {
    throw std::runtime_error(
        "SeedFinderOptions not in ACTS internal units in BinnedSPGroup");
  }
  static_assert(
      std::is_same<typename container_t::value_type,
                   const external_spacepoint_t*>::value,
      "Container does not contain type this class was templated with");

  // get region of interest (or full detector if configured accordingly)
  float phiMin = config.phiMin;
  float phiMax = config.phiMax;
  float zMin = config.zMin;
  float zMax = config.zMax;

  auto start = std::chrono::steady_clock::now();

  // sort by radius
  // add magnitude of beamPos to rMax to avoid excluding measurements
  // create number of bins equal to number of millimeters rMax
  // (worst case minR: configured minR + 1mm)
  // binSizeR allows to increase or reduce numRBins if needed
  size_t numRBins = static_cast<size_t>((config.rMax + options.beamPos.norm()) /
                                        config.binSizeR);

  // initialize original SP index locations
  std::vector<std::size_t> indexes(spacePoints.size());
  for (std::size_t i(0); i < indexes.size(); ++i) {
    indexes[i] = i;
  }

  // sort the SP indexes
  std::sort(indexes.begin(), indexes.end(),
            [&spacePoints](const std::size_t& a, const std::size_t& b) {
              return spacePoints[a]->r() < spacePoints[b]->r();
            });

  // fill the grid with space points sorted in radius
  // automatically all sps in all (z, phi) bins will be sorted in R
  for (const auto& idx : indexes) {
    const external_spacepoint_t& sp = *spacePoints[idx];

    const auto& [spPosition, variance] =
        toGlobal(sp, config.zAlign, config.rAlign, config.sigmaError);

    float spX = spPosition[0];
    float spY = spPosition[1];
    float spZ = spPosition[2];

    // store x,y,z values in extent
    rRangeSPExtent.extend({spX, spY, spZ});

    if (spZ > zMax || spZ < zMin) {
      continue;
    }
    float spPhi = std::atan2(spY, spX);
    if (spPhi > phiMax || spPhi < phiMin) {
      continue;
    }

    auto isp = std::make_unique<InternalSpacePoint<external_spacepoint_t>>(
        sp, spPosition, options.beamPos, variance);
    // calculate r-Bin index and protect against overflow (underflow not
    // possible)
    size_t rIndex = static_cast<size_t>(isp->radius() / config.binSizeR);
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins) {
      continue;
    }

    // fill rbins into grid
    Acts::Vector2 spLocation(isp->phi(), isp->z());
    std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>&
        rbin = grid->atPosition(spLocation);
    rbin.push_back(std::move(isp));
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::chrono::milliseconds elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_seconds);
  std::cout << "elapsed time: " << elapsed_seconds.count() << " micro sec \n";

  m_binnedSP = std::move(grid);
  m_bottomBinFinder = botBinFinder;
  m_topBinFinder = tBinFinder;

  m_bins = config.zBinsCustomLooping;
}
