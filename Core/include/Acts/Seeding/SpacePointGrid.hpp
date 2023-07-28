// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <memory>

namespace Acts {

struct SpacePointGridConfig {
  // minimum pT to be found by seedFinder
  float minPt = 0;
  // maximum extension of sensitive detector layer relevant for seeding as
  // distance from x=y=0 (i.e. in r)
  float rMax = 0;
  // maximum extension of sensitive detector layer relevant for seeding in
  // positive direction in z
  float zMax = 0;
  // maximum extension of sensitive detector layer relevant for seeding in
  // negative direction in z
  float zMin = 0;
  // maximum distance in r from middle space point to bottom or top spacepoint
  float deltaRMax = 0;
  // maximum forward direction expressed as cot(theta)
  float cotThetaMax = 0;
  // maximum impact parameter in mm
  float impactMax = 0;
  // minimum phi value for phiAxis construction
  float phiMin = -M_PI;
  // maximum phi value for phiAxis construction
  float phiMax = M_PI;
  // Multiplicator for the number of phi-bins. The minimum number of phi-bins
  // depends on min_pt, magnetic field: 2*M_PI/(minPT particle phi-deflection).
  // phiBinDeflectionCoverage is a multiplier for this number. If
  // numPhiNeighbors (in the configuration of the BinFinders) is configured to
  // return 1 neighbor on either side of the current phi-bin (and you want to
  // cover the full phi-range of minPT), leave this at 1.
  int phiBinDeflectionCoverage = 1;
  // enable non equidistant binning in z
  std::vector<float> zBinEdges;
  bool isInInternalUnits = false;
  SpacePointGridConfig toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SpacePointGridConfig");
    }
    using namespace Acts::UnitLiterals;
    SpacePointGridConfig config = *this;
    config.isInInternalUnits = true;
    config.minPt /= 1_MeV;
    config.rMax /= 1_mm;
    config.zMax /= 1_mm;
    config.zMin /= 1_mm;
    config.deltaRMax /= 1_mm;

    return config;
  }
};

struct SpacePointGridOptions {
  // magnetic field
  float bFieldInZ = 0;
  bool isInInternalUnits = false;
  SpacePointGridOptions toInternalUnits() const {
    if (isInInternalUnits) {
      throw std::runtime_error(
          "Repeated conversion to internal units for SpacePointGridOptions");
    }
    using namespace Acts::UnitLiterals;
    SpacePointGridOptions options = *this;
    options.isInInternalUnits = true;
    options.bFieldInZ /= 1000_T;

    return options;
  }
};

template <typename external_spacepoint_t>
using SpacePointGrid = detail::Grid<
    std::vector<std::unique_ptr<InternalSpacePoint<external_spacepoint_t>>>,
    detail::Axis<detail::AxisType::Equidistant,
                 detail::AxisBoundaryType::Closed>,
    detail::Axis<detail::AxisType::Variable, detail::AxisBoundaryType::Bound>>;

class SpacePointGridCreator {
 public:
  template <typename external_spacepoint_t>
  static std::unique_ptr<SpacePointGrid<external_spacepoint_t>> createGrid(
      const Acts::SpacePointGridConfig& _config,
      const Acts::SpacePointGridOptions& _options);
};
}  // namespace Acts
#include "Acts/Seeding/SpacePointGrid.ipp"
