// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFinderConfigNA60.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <vector>

namespace Acts {

/// Planar Space Point bin is a 2D grid with (x, z) bins
/// It stores a vector of internal space points to external space points
template <typename external_spacepoint_t>
using PlanarSpacePointGrid = Acts::Grid<
    std::vector<
        std::unique_ptr<Acts::InternalSpacePoint<external_spacepoint_t>>>,
    Acts::Axis<Acts::AxisType::Variable, Acts::AxisBoundaryType::Bound>,
    Acts::Axis<Acts::AxisType::Variable, Acts::AxisBoundaryType::Bound>>;

/// Planar Binned Group
template <typename external_spacepoint_t>
using PlanarBinnedGroup =
    Acts::BinnedGroup<Acts::PlanarSpacePointGrid<external_spacepoint_t>>;

template <typename external_spacepoint_t>
using PlanarBinnedGroupIterator = Acts::BinnedGroupIterator<
    Acts::PlanarSpacePointGrid<external_spacepoint_t>>;

struct PlanarSpacePointGridConfig {
  // maximum extension of sensitive detector layer relevant for seeding in
  // positive direction in x
  float xMax = 0;
  // maximum extension of sensitive detector layer relevant for seeding in
  // negative direction in x
  float xMin = 0;
  // maximum extension of sensitive detector layer relevant for seeding in
  // positive direction in z
  float zMax = 0;
  // maximum extension of sensitive detector layer relevant for seeding in
  // negative direction in z
  float zMin = 0;
  // enable non equidistant binning in z
  std::vector<float> zBinEdges;
  // enable non equidistant binning in x
  std::vector<float> xBinEdges;
  bool isInInternalUnits = false;
  PlanarSpacePointGridConfig toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for "
          "PlanarSpacePointGridConfig");
    }
    using namespace Acts::UnitLiterals;
    PlanarSpacePointGridConfig config = *this;
    config.isInInternalUnits = true;
    config.xMax /= 1_mm;
    config.xMin /= 1_mm;
    config.zMax /= 1_mm;
    config.zMin /= 1_mm;

    return config;
  }

};

struct PlanarSpacePointGridOptions {
  // magnetic field
  float bFieldInZ = 0;
  bool isInInternalUnits = false;
  PlanarSpacePointGridOptions toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for "
          "PlanarSpacePointGridOptions");
    }
    using namespace Acts::UnitLiterals;
    PlanarSpacePointGridOptions options = *this;
    options.isInInternalUnits = true;
    options.bFieldInZ /= 1000_T;

    return options;
    }
};
  

/// Instructions on how to create and fill this grid specialization
class PlanarSpacePointGridCreator {
 public:
  template <typename external_spacepoint_t>
  static Acts::PlanarSpacePointGrid<external_spacepoint_t> createGrid(
      const Acts::PlanarSpacePointGridConfig& _config,
      const Acts::PlanarSpacePointGridOptions& _options);

  template <typename external_spacepoint_t,
            typename external_spacepoint_iterator_t, typename callable_t>
  static void fillGrid(
      const Acts::SeedFinderConfigNA60<external_spacepoint_t>& config,
      const Acts::SeedFinderOptionsNA60& options,
      Acts::PlanarSpacePointGrid<external_spacepoint_t>& grid,
      external_spacepoint_iterator_t spBegin,
      external_spacepoint_iterator_t spEnd, callable_t&& toGlobal,
      Acts::Extent& rRangeSPExtent);
};

}  // namespace Acts

#include "Acts/Seeding/detail/PlanarSpacePointGrid.ipp"
