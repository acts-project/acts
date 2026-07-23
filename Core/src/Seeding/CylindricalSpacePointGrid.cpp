// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CylindricalSpacePointGrid.hpp"

#include "Acts/Seeding/detail/SpacePointGridPhiBinning.hpp"

namespace Acts {

CylindricalSpacePointGrid::CylindricalSpacePointGrid(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : m_cfg(config), m_logger(std::move(_logger)) {
  if (m_cfg.phiMin < -std::numbers::pi_v<float> ||
      m_cfg.phiMax > std::numbers::pi_v<float>) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid: phiMin (" + std::to_string(m_cfg.phiMin) +
        ") and/or phiMax (" + std::to_string(m_cfg.phiMax) +
        ") are outside the allowed phi range, defined as "
        "[-std::numbers::pi_v<float>, std::numbers::pi_v<float>]");
  }
  if (m_cfg.phiMin > m_cfg.phiMax) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid: phiMin is bigger then phiMax");
  }
  if (m_cfg.rMin > m_cfg.rMax) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid: rMin is bigger then rMax");
  }
  if (m_cfg.zMin > m_cfg.zMax) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid: zMin is bigger than zMax");
  }

  const int phiBins = detail::computeSpacePointGridPhiBins(
      {.minPt = m_cfg.minPt,
       .bFieldInZ = m_cfg.bFieldInZ,
       .rMax = m_cfg.rMax,
       .deltaRMax = m_cfg.deltaRMax,
       .impactMax = m_cfg.impactMax,
       .phiBinDeflectionCoverage = m_cfg.phiBinDeflectionCoverage,
       .maxPhiBins = m_cfg.maxPhiBins},
      logger());

  PhiAxisType phiAxis(AxisClosed, m_cfg.phiMin, m_cfg.phiMax, phiBins);

  // vector that will store the edges of the bins of z
  std::vector<double> zValues;

  // If zBinEdges is not defined, calculate the edges as zMin + bin * zBinSize
  if (m_cfg.zBinEdges.empty()) {
    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering
    const float zBinSize = m_cfg.cotThetaMax * m_cfg.deltaRMax;
    const float zBins =
        std::max(1.f, std::floor((m_cfg.zMax - m_cfg.zMin) / zBinSize));

    zValues.reserve(static_cast<int>(zBins));
    for (int bin = 0; bin <= static_cast<int>(zBins); bin++) {
      const double edge =
          m_cfg.zMin + bin * ((m_cfg.zMax - m_cfg.zMin) / zBins);
      zValues.push_back(edge);
    }
  } else {
    // Use the zBinEdges defined in the m_cfg
    zValues.reserve(m_cfg.zBinEdges.size());
    for (float bin : m_cfg.zBinEdges) {
      zValues.push_back(bin);
    }
  }

  std::vector<double> rValues;
  rValues.reserve(std::max(2ul, m_cfg.rBinEdges.size()));
  if (m_cfg.rBinEdges.empty()) {
    rValues = {m_cfg.rMin, m_cfg.rMax};
  } else {
    rValues.insert(rValues.end(), m_cfg.rBinEdges.begin(),
                   m_cfg.rBinEdges.end());
  }

  ZAxisType zAxis(AxisOpen, std::move(zValues));
  RAxisType rAxis(AxisOpen, std::move(rValues));

  ACTS_VERBOSE("Defining Grid:");
  ACTS_VERBOSE("- Phi Axis: " << phiAxis);
  ACTS_VERBOSE("- Z axis  : " << zAxis);
  ACTS_VERBOSE("- R axis  : " << rAxis);

  GridType grid(
      std::make_tuple(std::move(phiAxis), std::move(zAxis), std::move(rAxis)));
  m_binnedGroup.emplace(std::move(grid), m_cfg.bottomBinFinder.value(),
                        m_cfg.topBinFinder.value(), m_cfg.navigation);
  m_grid = &m_binnedGroup->grid();
}

void CylindricalSpacePointGrid::clear() {
  for (std::size_t i = 0; i < grid().size(); ++i) {
    BinType& bin = grid().at(i);
    bin.clear();
  }
  m_counter = 0;
}

std::optional<std::size_t> CylindricalSpacePointGrid::insert(
    SpacePointIndex index, float phi, float z, float r) {
  const std::optional<std::size_t> gridIndex = binIndex(phi, z, r);
  if (gridIndex.has_value()) {
    BinType& bin = grid().at(*gridIndex);
    bin.push_back(index);
    ++m_counter;
  }
  return gridIndex;
}

void CylindricalSpacePointGrid::extend(
    const SpacePointContainer::ConstRange& spacePoints) {
  ACTS_VERBOSE("Inserting " << spacePoints.size()
                            << " space points to the grid");

  for (const ConstSpacePointProxy& sp : spacePoints) {
    insert(sp);
  }
}

void CylindricalSpacePointGrid::sortBinsByR(
    const SpacePointContainer& spacePoints) {
  ACTS_VERBOSE("Sorting the grid");

  for (std::size_t i = 0; i < grid().size(); ++i) {
    BinType& bin = grid().at(i);
    std::ranges::sort(bin, {}, [&](SpacePointIndex spIndex) {
      return spacePoints[spIndex].zr()[1];
    });
  }

  ACTS_VERBOSE(
      "Number of space points inserted (within grid range): " << m_counter);
}

Range1D<float> CylindricalSpacePointGrid::computeRadiusRange(
    const SpacePointContainer& spacePoints) const {
  float minRange = std::numeric_limits<float>::max();
  float maxRange = std::numeric_limits<float>::lowest();
  for (const BinType& bin : grid()) {
    if (bin.empty()) {
      continue;
    }
    auto first = spacePoints[bin.front()];
    auto last = spacePoints[bin.back()];
    minRange = std::min(first.zr()[1], minRange);
    maxRange = std::max(last.zr()[1], maxRange);
  }
  return {minRange, maxRange};
}

}  // namespace Acts
