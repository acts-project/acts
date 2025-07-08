// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/CylindricalSpacePointGrid2.hpp"

namespace Acts::Experimental {

CylindricalSpacePointGrid2::CylindricalSpacePointGrid2(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : m_cfg(config), m_logger(std::move(_logger)) {
  if (m_cfg.phiMin < -std::numbers::pi_v<float> ||
      m_cfg.phiMax > std::numbers::pi_v<float>) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid2: phiMin (" + std::to_string(m_cfg.phiMin) +
        ") and/or phiMax (" + std::to_string(m_cfg.phiMax) +
        ") are outside the allowed phi range, defined as "
        "[-std::numbers::pi_v<float>, std::numbers::pi_v<float>]");
  }
  if (m_cfg.phiMin > m_cfg.phiMax) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid2: phiMin is bigger then phiMax");
  }
  if (m_cfg.rMin > m_cfg.rMax) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid2: rMin is bigger then rMax");
  }
  if (m_cfg.zMin > m_cfg.zMax) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid2: zMin is bigger than zMax");
  }

  int phiBins = 0;
  if (m_cfg.bFieldInZ == 0) {
    // for no magnetic field, use the maximum number of phi bins
    phiBins = m_cfg.maxPhiBins;
    ACTS_VERBOSE(
        "B-Field is 0 (z-coordinate), setting the number of bins in phi to "
        << phiBins);
  } else {
    // calculate circle intersections of helix and max detector radius in mm.
    // bFieldInZ is in (pT/radius) natively, no need for conversion
    float minHelixRadius = m_cfg.minPt / m_cfg.bFieldInZ;

    // sanity check: if yOuter takes the square root of a negative number
    if (minHelixRadius < m_cfg.rMax * 0.5) {
      throw std::domain_error(
          "The value of minHelixRadius cannot be smaller than rMax / 2. Please "
          "check the m_cfguration of bFieldInZ and minPt");
    }

    float maxR2 = m_cfg.rMax * m_cfg.rMax;
    float xOuter = maxR2 / (2 * minHelixRadius);
    float yOuter = std::sqrt(maxR2 - xOuter * xOuter);
    float outerAngle = std::atan(xOuter / yOuter);
    // intersection of helix and max detector radius minus maximum R distance
    // from middle SP to top SP
    float innerAngle = 0;
    float rMin = m_cfg.rMax;
    if (m_cfg.rMax > m_cfg.deltaRMax) {
      rMin = m_cfg.rMax - m_cfg.deltaRMax;
      float innerCircleR2 =
          (m_cfg.rMax - m_cfg.deltaRMax) * (m_cfg.rMax - m_cfg.deltaRMax);
      float xInner = innerCircleR2 / (2 * minHelixRadius);
      float yInner = std::sqrt(innerCircleR2 - xInner * xInner);
      innerAngle = std::atan(xInner / yInner);
    }

    // evaluating the azimutal deflection including the maximum impact parameter
    float deltaAngleWithMaxD0 =
        std::abs(std::asin(m_cfg.impactMax / rMin) -
                 std::asin(m_cfg.impactMax / m_cfg.rMax));

    // evaluating delta Phi based on the inner and outer angle, and the azimutal
    // deflection including the maximum impact parameter
    // Divide by m_cfg.phiBinDeflectionCoverage since we combine
    // m_cfg.phiBinDeflectionCoverage number of consecutive phi bins in the
    // seed making step. So each individual bin should cover
    // 1/m_cfg.phiBinDeflectionCoverage of the maximum expected azimutal
    // deflection
    float deltaPhi = (outerAngle - innerAngle + deltaAngleWithMaxD0) /
                     m_cfg.phiBinDeflectionCoverage;

    // sanity check: if the delta phi is equal to or less than zero, we'll be
    // creating an infinite or a negative number of bins, which would be bad!
    if (deltaPhi <= 0.f) {
      throw std::domain_error(
          "Delta phi value is equal to or less than zero, leading to an "
          "impossible number of bins (negative or infinite)");
    }

    // divide 2pi by angle delta to get number of phi-bins
    // size is always 2pi even for regions of interest
    phiBins = static_cast<int>(std::ceil(2 * std::numbers::pi / deltaPhi));
    // need to scale the number of phi bins accordingly to the number of
    // consecutive phi bins in the seed making step.
    // Each individual bin should be approximately a fraction (depending on this
    // number) of the maximum expected azimutal deflection.

    // set protection for large number of bins, by default it is large
    phiBins = std::min(phiBins, m_cfg.maxPhiBins);
  }

  Axis phiAxis(AxisClosed, m_cfg.phiMin, m_cfg.phiMax, phiBins);

  // vector that will store the edges of the bins of z
  std::vector<double> zValues{};

  // If zBinEdges is not defined, calculate the edges as zMin + bin * zBinSize
  if (m_cfg.zBinEdges.empty()) {
    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering
    float zBinSize = m_cfg.cotThetaMax * m_cfg.deltaRMax;
    float zBins =
        std::max(1.f, std::floor((m_cfg.zMax - m_cfg.zMin) / zBinSize));

    zValues.reserve(static_cast<int>(zBins));
    for (int bin = 0; bin <= static_cast<int>(zBins); bin++) {
      double edge = m_cfg.zMin + bin * ((m_cfg.zMax - m_cfg.zMin) / zBins);
      zValues.push_back(edge);
    }
  } else {
    // Use the zBinEdges defined in the m_cfg
    zValues.reserve(m_cfg.zBinEdges.size());
    for (float bin : m_cfg.zBinEdges) {
      zValues.push_back(bin);
    }
  }

  std::vector<double> rValues{};
  rValues.reserve(std::max(2ul, m_cfg.rBinEdges.size()));
  if (m_cfg.rBinEdges.empty()) {
    rValues = {m_cfg.rMin, m_cfg.rMax};
  } else {
    rValues.insert(rValues.end(), m_cfg.rBinEdges.begin(),
                   m_cfg.rBinEdges.end());
  }

  Axis zAxis(AxisOpen, std::move(zValues));
  Axis rAxis(AxisOpen, std::move(rValues));

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

void CylindricalSpacePointGrid2::insert(const ConstSpacePointProxy2& sp) {
  Vector3 position(sp.phi(), sp.z(), sp.r());
  if (!grid().isInside(position)) {
    return;
  }

  std::size_t globIndex = grid().globalBinFromPosition(position);
  auto& bin = grid().at(globIndex);
  bin.push_back(sp.index());
  ++m_counter;
}

void CylindricalSpacePointGrid2::extend(
    const SpacePointContainer2::ConstRange& spacePoints) {
  ACTS_VERBOSE("Inserting " << spacePoints.size()
                            << " space points to the grid");

  for (const auto& sp : spacePoints) {
    insert(sp);
  }
}

void CylindricalSpacePointGrid2::fill(const SpacePointContainer2& spacePoints) {
  extend(spacePoints.range({0, spacePoints.size()}));
  sort(spacePoints);
}

void CylindricalSpacePointGrid2::sort(const SpacePointContainer2& spacePoints) {
  ACTS_VERBOSE("Sorting the grid");

  for (std::size_t i = 0; i < grid().size(); ++i) {
    auto& bin = grid().at(i);
    std::ranges::sort(bin, {}, [&](SpacePointIndex2 spIndex) {
      return spacePoints[spIndex].r();
    });
  }

  ACTS_VERBOSE(
      "Number of space points inserted (within grid range): " << m_counter);
}

Range1D<float> CylindricalSpacePointGrid2::computeRadiusRange(
    const SpacePointContainer2& spacePoints) const {
  float minRange = std::numeric_limits<float>::max();
  float maxRange = std::numeric_limits<float>::lowest();
  for (const auto& coll : grid()) {
    if (coll.empty()) {
      continue;
    }
    auto first = spacePoints[coll.front()];
    auto last = spacePoints[coll.back()];
    minRange = std::min(first.r(), minRange);
    maxRange = std::max(last.r(), maxRange);
  }
  return {minRange, maxRange};
}

}  // namespace Acts::Experimental
