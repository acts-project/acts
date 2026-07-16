// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SphericalSpacePointGrid.hpp"

#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>

namespace Acts::Experimental {

SphericalSpacePointGrid::SphericalSpacePointGrid(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : m_cfg(config), m_logger(std::move(_logger)) {
  if (m_cfg.phiMin < -std::numbers::pi_v<float> ||
      m_cfg.phiMax > std::numbers::pi_v<float>) {
    throw std::runtime_error(
        "SphericalSpacePointGrid: phiMin (" + std::to_string(m_cfg.phiMin) +
        ") and/or phiMax (" + std::to_string(m_cfg.phiMax) +
        ") are outside the allowed phi range, defined as "
        "[-std::numbers::pi_v<float>, std::numbers::pi_v<float>]");
  }
  if (m_cfg.phiMin > m_cfg.phiMax) {
    throw std::runtime_error(
        "SphericalSpacePointGrid: phiMin is bigger then phiMax");
  }
  if (m_cfg.rMin > m_cfg.rMax) {
    throw std::runtime_error(
        "SphericalSpacePointGrid: rMin is bigger then rMax");
  }
  if (m_cfg.etaMin > m_cfg.etaMax) {
    throw std::runtime_error(
        "SphericalSpacePointGrid: etaMin is bigger than etaMax");
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
    const float minHelixRadius = m_cfg.minPt / m_cfg.bFieldInZ;

    // sanity check: if yOuter takes the square root of a negative number
    if (minHelixRadius < m_cfg.rMax * 0.5) {
      throw std::domain_error(
          "The value of minHelixRadius cannot be smaller than rMax / 2. Please "
          "check the configuration of bFieldInZ and minPt");
    }

    // x = rMax^2 / (2 * minHelixRadius)
    // y = cathetus(rMax, x)
    // R / x = 2 * minHelixRadius / R
    // outerAngle = x / y = x / sqrt(rMax^2 - x^2) = 1 / cath(rMax / x, 1)
    //            = 1 / cath(2 * minHelixRadius / rMax, 1)
    const float outerAngle =
        std::atan(1.f / fastCathetus(2 * minHelixRadius / m_cfg.rMax, 1));
    // intersection of helix and max detector radius minus maximum R distance
    // from middle SP to top SP
    float innerAngle = 0;
    float rMin = m_cfg.rMax;
    if (m_cfg.rMax > m_cfg.deltaRMax) {
      const float innerCircleR = m_cfg.rMax - m_cfg.deltaRMax;
      rMin = innerCircleR;
      innerAngle =
          std::atan(1.f / fastCathetus(2 * minHelixRadius / innerCircleR, 1));
    }

    // evaluating the azimutal deflection including the maximum impact
    // parameter. A track with |d0| >= r cannot reach radius r, so the azimuthal
    // deflection saturates at pi/2: clamp the asin arguments to [0, 1] so that
    // impactMax >= (rMax - deltaRMax) falls back to coarse / full-2pi phi
    // coverage instead of producing a NaN (and hence a zero phi-bin count and a
    // "Invalid binning" exception).
    const float sinInner = std::min(1.f, m_cfg.impactMax / rMin);
    const float sinOuter = std::min(1.f, m_cfg.impactMax / m_cfg.rMax);
    const float deltaAngleWithMaxD0 =
        std::abs(std::asin(sinInner) - std::asin(sinOuter));

    // evaluating delta Phi based on the inner and outer angle, and the azimutal
    // deflection including the maximum impact parameter
    // Divide by m_cfg.phiBinDeflectionCoverage since we combine
    // m_cfg.phiBinDeflectionCoverage number of consecutive phi bins in the
    // seed making step. So each individual bin should cover
    // 1/m_cfg.phiBinDeflectionCoverage of the maximum expected azimutal
    // deflection
    const float deltaPhi = (outerAngle - innerAngle + deltaAngleWithMaxD0) /
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

  PhiAxisType phiAxis(AxisClosed, m_cfg.phiMin, m_cfg.phiMax, phiBins);

  // vector that will store the edges of the eta bins. The edges are stored as
  // cot(theta) = sinh(eta), so a space point
  // bins by z / r into exactly the eta bin it belongs to.
  std::vector<double> etaValues;

  // If etaBinEdges is not defined, calculate equidistant eta edges.
  if (m_cfg.etaBinEdges.empty()) {
    const float etaBinSize = m_cfg.deltaEtaMax;
    const float etaBins =
        std::max(1.f, std::floor((m_cfg.etaMax - m_cfg.etaMin) / etaBinSize));

    etaValues.reserve(static_cast<int>(etaBins) + 1);
    for (int bin = 0; bin <= static_cast<int>(etaBins); bin++) {
      const double eta =
          m_cfg.etaMin + bin * ((m_cfg.etaMax - m_cfg.etaMin) / etaBins);
      etaValues.push_back(std::sinh(eta));
    }
  } else {
    // Use the etaBinEdges defined in the m_cfg, mapped to cot(theta).
    etaValues.reserve(m_cfg.etaBinEdges.size());
    for (float eta : m_cfg.etaBinEdges) {
      etaValues.push_back(std::sinh(eta));
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

  CotThetaAxisType cotThetaAxis(AxisOpen, std::move(etaValues));
  RAxisType rAxis(AxisOpen, std::move(rValues));

  ACTS_VERBOSE("Defining Grid:");
  ACTS_VERBOSE("- Phi Axis: " << phiAxis);
  ACTS_VERBOSE("- cot(theta) axis (sinh(eta) edges): " << cotThetaAxis);
  ACTS_VERBOSE("- R axis  : " << rAxis);

  GridType grid(std::make_tuple(std::move(phiAxis), std::move(cotThetaAxis),
                                std::move(rAxis)));
  m_binnedGroup.emplace(std::move(grid), m_cfg.bottomBinFinder.value(),
                        m_cfg.topBinFinder.value(), m_cfg.navigation);
  m_grid = &m_binnedGroup->grid();
}

void SphericalSpacePointGrid::clear() {
  for (std::size_t i = 0; i < grid().size(); ++i) {
    BinType& bin = grid().at(i);
    bin.clear();
  }
  m_counter = 0;
}

std::optional<std::size_t> SphericalSpacePointGrid::insert(
    SpacePointIndex index, float phi, float cotTheta, float r) {
  const std::optional<std::size_t> gridIndex = binIndex(phi, cotTheta, r);
  if (gridIndex.has_value()) {
    BinType& bin = grid().at(*gridIndex);
    bin.push_back(index);
    ++m_counter;
  }
  return gridIndex;
}

void SphericalSpacePointGrid::extend(
    const SpacePointContainer::ConstRange& spacePoints) {
  ACTS_VERBOSE("Inserting " << spacePoints.size()
                            << " space points to the grid");

  for (const ConstSpacePointProxy& sp : spacePoints) {
    insert(sp);
  }
}

void SphericalSpacePointGrid::sortBinsByR(
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

Range1D<float> SphericalSpacePointGrid::computeRadiusRange(
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

}  // namespace Acts::Experimental
