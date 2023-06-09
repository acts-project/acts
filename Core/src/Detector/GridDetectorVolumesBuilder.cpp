// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/GridDetectorVolumesBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/GridDetectorVolumesBuilder.hpp"
#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <set>
#include <stdexcept>

Acts::Experimental::GridDetectorVolumesBuilder::GridDetectorVolumesBuilder(
    const Acts::Experimental::GridDetectorVolumesBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> mlogger)
    : IDetectorComponentBuilder(), m_cfg(cfg), m_logger(std::move(mlogger)) {
  // Early bail out
  if (m_cfg.binning.size() != 3u) {
    throw std::invalid_argument(
        "GridDetectorVolumesBuilder: Invalid binning, exact 3 binning values "
        "need to be provided, either z/r/phi or x/y/z.");
  }

  // Sort & check the configuration
  std::sort(m_cfg.binning.begin(), m_cfg.binning.end(),
            [&](const ProtoBinning& a, const ProtoBinning& b) {
              return a.binValue < b.binValue;
            });

  m_binningValues = {m_cfg.binning[0].binValue, m_cfg.binning[1].binValue,
                     m_cfg.binning[2].binValue};

  // Unknown binning
  if (m_binningValues != cylindricalBinning and
      m_binningValues != cartesianBinning) {
    throw std::invalid_argument(
        "GridDetectorVolumesBuilder: Invalid binning, must be either z/r/phi "
        "or x/y/z");
  }
  // Closed disallowed and catch numerical phi closure
  std::size_t closedDisallowed =
      (m_binningValues == cartesianBinning) ? 3u : 2u;
  for (auto [ib, b] : enumerate(m_cfg.binning)) {
    if (ib < closedDisallowed and
        b.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
      throw std::invalid_argument(
          "GridDetectorVolumesBuilder: Invalid binning, only for z/r/phi the "
          "last binning can be closed.");
    }
    if (b.binValue == binPhi) {
      // Phi binning close to 2pi - renormalize
      if (std::abs(b.edges.back() - b.edges.front() - 2 * M_PI) <
          m_cfg.phiClosedTolerance) {
        m_closedPhiSetup = true;
        ACTS_DEBUG("Closed phi setup detected, adjusting binning to [-pi,pi]");
      }
    }
  }
}

Acts::Experimental::DetectorComponent
Acts::Experimental::GridDetectorVolumesBuilder::construct(
    const Acts::GeometryContext& gctx) const {
  // Screen output grid type
  std::string gridType =
      (m_binningValues == cylindricalBinning) ? "cylindrical" : "cartesian";

  ACTS_DEBUG("Constructing " << gridType << " grid setup " << m_cfg.name);

  std::array<Acts::detail::AxisType, 3u> axesType = {
      {m_cfg.binning[0].axisType, m_cfg.binning[1].axisType,
       m_cfg.binning[2].axisType}};

  // Get into the different cases: equidistant x 3
  if (axesType == eeeType and m_closedPhiSetup) {
    detail::GridAxisGenerators::EqBoundEqBoundEqClosed eBeBeC{
        {m_cfg.binning[0].edges.front(), m_cfg.binning[0].edges.back()},
        m_cfg.binning[0].bins(),
        {m_cfg.binning[1].edges.front(), m_cfg.binning[1].edges.back()},
        m_cfg.binning[1].bins(),
        {-M_PI, M_PI},
        m_cfg.binning[2].bins()};
    // With closed phi
    return constructT(gctx, eBeBeC);
  } else if (axesType == eeeType) {
    detail::GridAxisGenerators::EqBoundEqBoundEqBound eBeBeB{
        {m_cfg.binning[0].edges.front(), m_cfg.binning[0].edges.back()},
        m_cfg.binning[0].bins(),
        {m_cfg.binning[1].edges.front(), m_cfg.binning[1].edges.back()},
        m_cfg.binning[1].bins(),
        {m_cfg.binning[2].edges.front(), m_cfg.binning[2].edges.back()},
        m_cfg.binning[2].bins()};
    // With triple bound
    return constructT(gctx, eBeBeB);
  }

  // Get into the different cases: equidistant x 3
  if (axesType == vveType and m_closedPhiSetup) {
    detail::GridAxisGenerators::VarBoundVarBoundEqClosed vBvBeC{
        m_cfg.binning[0].edges,
        m_cfg.binning[1].edges,
        {-M_PI, M_PI},
        m_cfg.binning[2].bins()};

    return constructT(gctx, vBvBeC);
  } else if (axesType == vveType) {
    detail::GridAxisGenerators::VarBoundVarBoundEqBound vBvBeB{
        m_cfg.binning[0].edges,
        m_cfg.binning[1].edges,
        {m_cfg.binning[2].edges.front(), m_cfg.binning[2].edges.back()},
        m_cfg.binning[2].bins()};
    return constructT(gctx, vBvBeB);
  }

  // Get into the different cases: variable x 3
  if (axesType == vveType) {
    detail::GridAxisGenerators::VarBoundVarBoundVarBound vBvBvB{
        m_cfg.binning[0].edges, m_cfg.binning[1].edges, m_cfg.binning[2].edges};

    return constructT(gctx, vBvBvB);
  }

  return Acts::Experimental::DetectorComponent{
      {}, {}, RootDetectorVolumes{{}, tryNoVolumes()}};
}
