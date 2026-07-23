// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SphericalSpacePointGrid.hpp"

#include "Acts/Seeding/detail/SpacePointGridPhiBinning.hpp"

#include <cmath>

namespace Acts::Experimental {

namespace {

float radiusProjection(const ConstSpacePointProxy& sp) {
  return sp.zr()[1];
}

}  // namespace

SphericalSpacePointGrid::SphericalSpacePointGrid(
    const Config& config, std::unique_ptr<const Logger> _logger)
    : GridBase(std::move(_logger)), m_cfg(config) {
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

  const int phiBins = Acts::detail::computeSpacePointGridPhiBins(
      {.minPt = m_cfg.minPt,
       .bFieldInZ = m_cfg.bFieldInZ,
       .rMax = m_cfg.rMax,
       .deltaRMax = m_cfg.deltaRMax,
       .impactMax = m_cfg.impactMax,
       .phiBinDeflectionCoverage = m_cfg.phiBinDeflectionCoverage,
       .maxPhiBins = m_cfg.maxPhiBins},
      logger());

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
  initializeGrid(std::move(grid), m_cfg.bottomBinFinder.value(),
                 m_cfg.topBinFinder.value(), m_cfg.navigation);
}

void SphericalSpacePointGrid::sortBinsByR(
    const SpacePointContainer& spacePoints) {
  sortBinsBy(spacePoints, radiusProjection);
}

Range1D<float> SphericalSpacePointGrid::computeRadiusRange(
    const SpacePointContainer& spacePoints) const {
  return computeRange(spacePoints, radiusProjection);
}

}  // namespace Acts::Experimental
