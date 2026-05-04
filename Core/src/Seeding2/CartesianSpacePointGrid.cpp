// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/CartesianSpacePointGrid.hpp"

#include "Acts/Definitions/Common.hpp"

namespace Acts {

CartesianSpacePointGrid::CartesianSpacePointGrid(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : m_cfg(config), m_logger(std::move(_logger)) {
  if (m_cfg.xMin > m_cfg.xMax) {
    throw std::invalid_argument(
        "CartesianSpacePointGrid: xMin is bigger then xMax");
  }
  if (m_cfg.yMin > m_cfg.yMax) {
    throw std::invalid_argument(
        "CartesianSpacePointGrid: yMin is bigger then yMax");
  }
  if (m_cfg.zMin > m_cfg.zMax) {
    throw std::invalid_argument(
        "CartesianSpacePointGrid: zMin is bigger than zMax");
  }
  if (m_cfg.sortingCoord == Acts::eTime) {
    throw std::invalid_argument(
        "CartesianSpacePointGrid: time is not supported as a spatial sorting "
        "dimension");
  }

  XAxisType xAxis(m_cfg.xMin, m_cfg.xMax, m_cfg.nXbins);
  YAxisType yAxis(m_cfg.yMin, m_cfg.yMax, m_cfg.nYbins);
  YAxisType zAxis(m_cfg.zMin, m_cfg.zMax, m_cfg.nZbins);

  ACTS_VERBOSE("Defining Grid:");
  ACTS_VERBOSE("- x Axis  : " << xAxis);
  ACTS_VERBOSE("- y Axis  : " << yAxis);
  ACTS_VERBOSE("- z Axis  : " << zAxis);

  GridType grid(
      std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)));
  m_binnedGroup.emplace(std::move(grid), m_cfg.bottomBinFinder.value(),
                        m_cfg.topBinFinder.value(), m_cfg.navigation);
  m_grid = &m_binnedGroup->grid();
  m_sortCoordGetter = [c = m_cfg.sortingCoord, d = m_cfg.sortingDirection](
                          const ConstSpacePointProxy2& p) {
    return d * p.xyz()[c];
  };
}

void CartesianSpacePointGrid::clear() {
  for (std::size_t i = 0; i < grid().size(); ++i) {
    BinType& bin = grid().at(i);
    bin.clear();
  }
  m_counter = 0;
}

std::optional<std::size_t> CartesianSpacePointGrid::insert(
    SpacePointIndex index, float x, float y, float z) {
  const std::optional<std::size_t> gridIndex = binIndex(x, y, z);
  if (gridIndex.has_value()) {
    BinType& bin = grid().at(*gridIndex);
    bin.push_back(index);
    ++m_counter;
  }
  return gridIndex;
}

void CartesianSpacePointGrid::extend(
    const SpacePointContainer2::ConstRange& spacePoints) {
  ACTS_VERBOSE("Inserting " << spacePoints.size()
                            << " space points to the grid");

  for (const ConstSpacePointProxy2& sp : spacePoints) {
    insert(sp);
  }
}

void CartesianSpacePointGrid::sortBinsByCoord(
    const SpacePointContainer2& spacePoints) {
  ACTS_VERBOSE("Sorting the grid");

  for (std::size_t i = 0; i < grid().size(); ++i) {
    BinType& bin = grid().at(i);
    std::ranges::sort(bin, {}, [&](SpacePointIndex2 spIndex) {
      return m_sortCoordGetter(spacePoints[spIndex]);
    });
  }

  ACTS_VERBOSE(
      "Number of space points inserted (within grid range): " << m_counter);
}

Range1D<float> CartesianSpacePointGrid::computeCoordRange(
    const SpacePointContainer2& spacePoints) const {
  float minRange = std::numeric_limits<float>::max();
  float maxRange = std::numeric_limits<float>::lowest();
  for (const BinType& bin : grid()) {
    if (bin.empty()) {
      continue;
    }
    auto first = spacePoints[bin.front()];
    auto last = spacePoints[bin.back()];
    minRange = std::min(m_sortCoordGetter(first), minRange);
    maxRange = std::max(m_sortCoordGetter(last), maxRange);
  }
  return {minRange, maxRange};
}

}  // namespace Acts
