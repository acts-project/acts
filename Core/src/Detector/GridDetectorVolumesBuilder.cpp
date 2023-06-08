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
    std::unique_ptr<const Acts::Logger> logger)
    : IDetectorComponentBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  // Early bail out
  if (m_cfg.binning.size() != 3u) {
    throw std::invalid_argument(
        "GridDetectorVolumesBuilder: Invalid binning, exact 3 binning values "
        "need to "
        "be provided, either z/r/phi or x/y/z.");
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

  // Check that only in the (z,r,phi) binning the last value can be closed
  if (m_binningValues == cartesianBinning) {
    for (const auto& b : m_cfg.binning) {
      if (b.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
        throw std::invalid_argument(
            "GridDetectorVolumesBuilder: Invalid binning, only in z/r/phi the "
            "last bin can be closed.");
      }
    }
  } else {
    // Check that only in the (z,r,phi) binning the last value can be closed
    for (auto [ib, b] : enumerate(m_cfg.binning)) {
      if (ib < 2u and
          b.boundaryType == Acts::detail::AxisBoundaryType::Closed) {
        throw std::invalid_argument(
            "GridDetectorVolumesBuilder: Invalid binning, only in z/r/phi the "
            "last bin can be closed.");
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
  if (axesType == eeeType) {
    detail::GridAxisGenerators::EqBoundEqBoundEqClosed eBeBeC{
        {m_cfg.binning[0].edges.front(), m_cfg.binning[0].edges.back()},
        m_cfg.binning[0].bins(),
        {m_cfg.binning[1].edges.front(), m_cfg.binning[1].edges.back()},
        m_cfg.binning[1].bins(),
        {-M_PI, M_PI},
        m_cfg.binning[2].bins()};

    return constructT(gctx, eBeBeC);
  }

  // Get into the different cases: equidistant x 3
  if (axesType == vveType) {
    detail::GridAxisGenerators::VarBoundVarBoundEqClosed vBvBeC{
        m_cfg.binning[0].edges,
        m_cfg.binning[1].edges,
        {-M_PI, M_PI},
        m_cfg.binning[2].bins()};

    return constructT(gctx, vBvBeC);
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
