// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
template <typename spacepoint_iterator_t>
Acts::BinnedSPGroup::BinnedSPGroup(
    spacepoint_iterator_t spBegin,
    spacepoint_iterator_t spEnd,
    std::shared_ptr<Acts::BinFinder> botBinFinder,
    std::shared_ptr<Acts::BinFinder> tBinFinder,
    std::unique_ptr<SpacePointGrid> grid,
    const SeedfinderConfig& _config) {
  auto config = _config.toInternalUnits();
  static_assert(
      std::is_same<
          typename std::iterator_traits<spacepoint_iterator_t>::value_type,
          Acts::SpacePoint*>::value,
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
  size_t numRBins = (config.rMax + config.beamPos.norm()) / config.binSizeR;
  std::vector<
      std::vector<Acts::SpacePoint*>>
      rBins(numRBins);
  for (spacepoint_iterator_t it = spBegin; it != spEnd; it++) {
    if (*it == nullptr) {
      continue;
    }
    Acts::SpacePoint* sp = *it;

    float spZ = sp->z();
    if (spZ > zMax || spZ < zMin) {
      continue;
    }

    float spPhi = sp->phi();
    if (spPhi > phiMax || spPhi < phiMin) {
      continue;
    }

    // calculate r-Bin index and protect against overflow (underflow not
    // possible)
    size_t rIndex = sp->radius() / config.binSizeR;
    // if index out of bounds, the SP is outside the region of interest
    if (rIndex >= numRBins) {
      continue;
    }
    rBins[rIndex].push_back(sp);
  }

  // if requested, it is possible to force sorting in R for each (z, phi) grid
  // bin
  if (config.forceRadialSorting) {
    for (auto& rbin : rBins) {
      std::sort(
          rbin.begin(), rbin.end(),
          [](Acts::SpacePoint* a,
             Acts::SpacePoint* b) {
            return a->radius() < b->radius();
          });
    }
  }

  // fill rbins into grid such that each grid bin is sorted in r
  // space points with delta r < rbin size can be out of order is sorting is not
  // requested
  for (auto& rbin : rBins) {
    for (auto& sp : rbin) {
      Acts::Vector2 spLocation(sp->phi(), sp->z());
      std::vector<Acts::SpacePoint*>&
          bin = grid->atPosition(spLocation);
      bin.push_back(std::move(sp));
    }
  }
  m_binnedSP = std::move(grid);
  m_bottomBinFinder = botBinFinder;
  m_topBinFinder = tBinFinder;

  m_bins = config.zBinsCustomLooping;
}
