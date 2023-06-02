// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/GridDetectorBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/GridDetectorBuilder.hpp"
#include "Acts/Detector/detail/CylindricalGridVolumesHelper.hpp"
#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <set>
#include <stdexcept>

Acts::Experimental::GridDetectorBuilder::GridDetectorBuilder(
    const Acts::Experimental::GridDetectorBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IDetectorBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  // Early bail out
  if (m_cfg.binning.size() != 3u) {
    throw std::invalid_argument(
        "GridDetectorBuilder: Invalid binning, exact 3 binning values need to "
        "be provided, either z/r/phi or x/y/z.");
  }

  // Sort & check the configuration
  std::sort(m_cfg.binning.begin(), m_cfg.binning.end(),
            [&](const Binning& a, const Binning& b) {
              return a.data.binvalue < b.data.binvalue;
            });

  m_binningValues = {m_cfg.binning[0].data.binvalue,
                     m_cfg.binning[1].data.binvalue,
                     m_cfg.binning[2].data.binvalue};
  // Unknown binning
  if (m_binningValues != cylindricalBinning and
      m_binningValues != cartesianBinning) {
    throw std::invalid_argument(
        "GridDetectorBuilder: Invalid binning, must be either z/r/phi or "
        "x/y/z");
  }

  // Check that only in the (z,r,phi) binning the last value can be closed
  if (m_binningValues == cartesianBinning) {
    for (const auto& b : m_cfg.binning) {
      if (b.data.option == Acts::BinningOption::closed) {
        throw std::invalid_argument(
            "GridDetectorBuilder: Invalid binning, only in z/r/phi the last "
            "bin can be closed.");
      }
    }
  } else {
    // Check that only in the (z,r,phi) binning the last value can be closed
    for (auto [ib, b] : enumerate(m_cfg.binning)) {
      if (ib < 2u and b.data.option == Acts::BinningOption::closed) {
        throw std::invalid_argument(
            "GridDetectorBuilder: Invalid binning, only in z/r/phi the last "
            "bin can be closed.");
      }
    }
  }
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::GridDetectorBuilder::construct(
    const Acts::GeometryContext& gctx) const {
  // Screen output grid type
  std::string gridType =
      (m_binningValues == cylindricalBinning) ? "cylindrical" : "cartesian";

  ACTS_DEBUG("Constructing " << gridType << " grid setup " << m_cfg.name);

  detail::GridAxisGenerators::EqBoundEqBoundEqClosed eqB2eqC{
      {m_cfg.binning[0].data.min, m_cfg.binning[0].data.max},
      m_cfg.binning[0].data.bins(),
      {m_cfg.binning[1].data.min, m_cfg.binning[1].data.max},
      m_cfg.binning[1].data.bins(),
      {-M_PI, M_PI},
      m_cfg.binning[2].data.bins()};

  using GridType =
      typename decltype(eqB2eqC)::template grid_type<std::vector<std::size_t>>;
  GridType grid(eqB2eqC());

  auto surfaces = m_cfg.surfaces();
  ACTS_DEBUG("Retreived " << surfaces.size()
                          << " surfaces to assign to the grid.");

  // We use the Indexed Surface grid functionality to assign surfaces to the
  IndexedSurfacesImpl<GridType> isGrid(std::move(grid), m_binningValues,
                                       m_cfg.transform);
  detail::IndexedGridFiller isGridFiller;
  detail::PolyhedronReferenceGenerator polyGen{true, m_cfg.polyhedronSegements};
  isGridFiller.fill(gctx, isGrid, surfaces, polyGen);

  using RootVolumesGridType =
      typename decltype(eqB2eqC)::template grid_type<std::size_t>;

  RootVolumesGridType rootVolumesGrid(eqB2eqC());

  auto gridVolumes =
      detail::CylindricalGridVolumesHelper::buildVolumes<GridType,
                                                         RootVolumesGridType>(
          gctx, isGrid.grid, rootVolumesGrid, surfaces, {});

  ACTS_DEBUG("Retreived " << gridVolumes.size() << " grid volumes.");

  // Create the detector
  return Detector::makeShared(
      m_cfg.name, gridVolumes,
      rootVolumesFromGrid(std::move(rootVolumesGrid), m_binningValues,
                          m_cfg.transform));
}
