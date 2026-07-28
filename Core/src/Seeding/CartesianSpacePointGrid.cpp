// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CartesianSpacePointGrid.hpp"

#include "Acts/Definitions/Common.hpp"

namespace Acts {

CartesianSpacePointGrid::CartesianSpacePointGrid(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : GridBase(std::move(_logger)), m_cfg(config) {
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
  ZAxisType zAxis(m_cfg.zMin, m_cfg.zMax, m_cfg.nZbins);

  ACTS_VERBOSE("Defining Grid:");
  ACTS_VERBOSE("- x Axis  : " << xAxis);
  ACTS_VERBOSE("- y Axis  : " << yAxis);
  ACTS_VERBOSE("- z Axis  : " << zAxis);

  GridType grid(
      std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)));
  initializeGrid(std::move(grid), m_cfg.bottomBinFinder.value(),
                 m_cfg.topBinFinder.value(), m_cfg.navigation);
  m_sortCoordGetter = [c = m_cfg.sortingCoord, d = m_cfg.sortingDirection](
                          const ConstSpacePointProxy& p) {
    return d * p.xyz()[c];
  };
}

void CartesianSpacePointGrid::sortBinsByCoord(
    const SpacePointContainer& spacePoints) {
  sortBinsBy(spacePoints, m_sortCoordGetter);
}

Range1D<float> CartesianSpacePointGrid::computeCoordRange(
    const SpacePointContainer& spacePoints) const {
  return computeRange(spacePoints, m_sortCoordGetter);
}

}  // namespace Acts
