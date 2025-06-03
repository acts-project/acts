// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/CylindricalSpacePointGrid2.hpp"

namespace Acts {

using namespace UnitLiterals;

CylindricalSpacePointGrid2::DerivedConfig
CylindricalSpacePointGrid2::Config::derive() const {
  DerivedConfig result;

  static_cast<Config&>(result) = *this;

  // TODO get rid of unit conversions
  {
    result.minPt /= 1_MeV;
    result.rMin /= 1_mm;
    result.rMax /= 1_mm;
    result.zMax /= 1_mm;
    result.zMin /= 1_mm;
    result.deltaRMax /= 1_mm;
    result.impactMax /= 1_mm;

    for (float& val : result.zBinEdges) {
      val /= 1_mm;
    }
    for (float& val : result.rBinEdges) {
      val /= 1_mm;
    }

    result.bFieldInZ /= 1000_T;
  }

  return result;
}

CylindricalSpacePointGrid2::CylindricalSpacePointGrid2(
    const DerivedConfig& config, std::unique_ptr<const Logger> _logger)
    : m_cfg(config), m_logger(std::move(_logger)) {
  if (m_cfg.phiMin < -std::numbers::pi_v<float> ||
      m_cfg.phiMax > std::numbers::pi_v<float>) {
    throw std::runtime_error(
        "CylindricalSpacePointGrid2: phiMin (" + std::to_string(m_cfg.phiMin) +
        ") and/or phiMax (" + std::to_string(m_cfg.phiMax) +
        ") are outside "
        "the allowed phi range, defined as "
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
  // for no magnetic field, create 100 phi-bins
  if (m_cfg.bFieldInZ == 0) {
    phiBins = 100;
    ACTS_VERBOSE(
        "B-Field is 0 (z-coordinate), setting the number of bins in phi to "
        << phiBins);
  } else {
    // calculate circle intersections of helix and max detector radius
    float minHelixRadius =
        config.minPt /
        (1_T * 1e6 *
         m_cfg.bFieldInZ);  // in mm -> R[mm] =pT[GeV] / (3·10−4×B[T])
                            // = pT[MeV] / (300 *Bz[kT])

    // sanity check: if yOuter takes the square root of a negative number
    if (minHelixRadius < config.rMax * 0.5) {
      throw std::domain_error(
          "The value of minHelixRadius cannot be smaller than rMax / 2. Please "
          "check the configuration of bFieldInZ and minPt");
    }

    float maxR2 = config.rMax * config.rMax;
    float xOuter = maxR2 / (2 * minHelixRadius);
    float yOuter = std::sqrt(maxR2 - xOuter * xOuter);
    float outerAngle = std::atan(xOuter / yOuter);
    // intersection of helix and max detector radius minus maximum R distance
    // from middle SP to top SP
    float innerAngle = 0;
    float rMin = config.rMax;
    if (config.rMax > config.deltaRMax) {
      rMin = config.rMax - config.deltaRMax;
      float innerCircleR2 =
          (config.rMax - config.deltaRMax) * (config.rMax - config.deltaRMax);
      float xInner = innerCircleR2 / (2 * minHelixRadius);
      float yInner = std::sqrt(innerCircleR2 - xInner * xInner);
      innerAngle = std::atan(xInner / yInner);
    }

    // evaluating the azimutal deflection including the maximum impact parameter
    float deltaAngleWithMaxD0 =
        std::abs(std::asin(config.impactMax / (rMin)) -
                 std::asin(config.impactMax / config.rMax));

    // evaluating delta Phi based on the inner and outer angle, and the azimutal
    // deflection including the maximum impact parameter
    // Divide by config.phiBinDeflectionCoverage since we combine
    // config.phiBinDeflectionCoverage number of consecutive phi bins in the
    // seed making step. So each individual bin should cover
    // 1/config.phiBinDeflectionCoverage of the maximum expected azimutal
    // deflection
    float deltaPhi = (outerAngle - innerAngle + deltaAngleWithMaxD0) /
                     config.phiBinDeflectionCoverage;

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
    phiBins = std::min(phiBins, config.maxPhiBins);
  }

  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> phiAxis(
      config.phiMin, config.phiMax, phiBins);

  // vector that will store the edges of the bins of z
  std::vector<double> zValues{};

  // If zBinEdges is not defined, calculate the edges as zMin + bin * zBinSize
  if (config.zBinEdges.empty()) {
    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering
    float zBinSize = config.cotThetaMax * config.deltaRMax;
    float zBins =
        std::max(1.f, std::floor((config.zMax - config.zMin) / zBinSize));

    zValues.reserve(static_cast<int>(zBins));
    for (int bin = 0; bin <= static_cast<int>(zBins); bin++) {
      double edge = config.zMin + bin * ((config.zMax - config.zMin) / zBins);
      zValues.push_back(edge);
    }

  } else {
    // Use the zBinEdges defined in the config
    zValues.reserve(config.zBinEdges.size());
    for (float bin : config.zBinEdges) {
      zValues.push_back(bin);
    }
  }

  std::vector<double> rValues{};
  rValues.reserve(std::max(2ul, config.rBinEdges.size()));
  if (config.rBinEdges.empty()) {
    rValues = {config.rMin, config.rMax};
  } else {
    rValues.insert(rValues.end(), config.rBinEdges.begin(),
                   config.rBinEdges.end());
  }

  Axis<AxisType::Variable, AxisBoundaryType::Open> zAxis(std::move(zValues));
  Axis<AxisType::Variable, AxisBoundaryType::Open> rAxis(std::move(rValues));

  ACTS_VERBOSE("Defining Grid:");
  ACTS_VERBOSE("- Phi Axis: " << phiAxis);
  ACTS_VERBOSE("- Z axis  : " << zAxis);
  ACTS_VERBOSE("- R axis  : " << rAxis);

  GridType grid(
      std::make_tuple(std::move(phiAxis), std::move(zAxis), std::move(rAxis)));
  m_binnedGroup.emplace(std::move(grid), m_cfg.bottomBinFinder.value(),
                        m_cfg.topBinFinder.value(), m_cfg.navigation);
  m_grid = &m_binnedGroup->grid();

  m_usedBinIndex = std::vector<bool>(m_grid->size(), false);
  m_rBinsIndex.reserve(m_grid->size());
}

void CylindricalSpacePointGrid2::extend(
    const SpacePointContainer2::ConstRange& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& phiColumn,
    const SpacePointContainer2::DenseColumn<float>& rColumn) {
  ACTS_VERBOSE("Inserting " << spacePoints.size()
                            << " space points to the grid");

  for (const auto& sp : spacePoints) {
    insert(sp, phiColumn, rColumn);
  }
}

void CylindricalSpacePointGrid2::fill(
    const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& phiColumn,
    const SpacePointContainer2::DenseColumn<float>& rColumn) {
  extend(spacePoints.range({0, spacePoints.size()}), phiColumn, rColumn);
  sort(spacePoints, rColumn);
}

void CylindricalSpacePointGrid2::sort(
    const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& rColumn) {
  ACTS_VERBOSE("Sorting the grid");

  for (std::size_t binIndex : m_rBinsIndex) {
    auto& rbin = grid().atPosition(binIndex);
    std::ranges::sort(rbin, {}, [&](const auto& sp) {
      return spacePoints.at(sp).extra(rColumn);
    });
  }

  m_usedBinIndex.assign(grid().size(), false);
  m_rBinsIndex.clear();

  ACTS_VERBOSE(
      "Number of space points inserted (within grid range): " << m_counter);
}

Range1D<float> CylindricalSpacePointGrid2::computeRadiusRange(
    const SpacePointContainer2& spacePoints,
    const SpacePointContainer2::DenseColumn<float>& rColumn) const {
  float minRange = std::numeric_limits<float>::max();
  float maxRange = std::numeric_limits<float>::lowest();
  for (const auto& coll : grid()) {
    if (coll.empty()) {
      continue;
    }
    auto first = spacePoints.at(coll.front());
    auto last = spacePoints.at(coll.back());
    minRange = std::min(first.extra(rColumn), minRange);
    maxRange = std::max(last.extra(rColumn), maxRange);
  }
  return {minRange, maxRange};
}

}  // namespace Acts
