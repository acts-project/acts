// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <memory>

namespace Acts {

struct SpacePointGridConfig {
  // magnetic field in kTesla
  float bFieldInZ;
  // minimum pT to be found by seedfinder in MeV
  float minPt;
  // maximum extension of sensitive detector layer relevant for seeding as
  // distance from x=y=0 (i.e. in r) in mm
  float rMax;
  // maximum extension of sensitive detector layer relevant for seeding in
  // positive direction in z in mm
  float zMax;
  // maximum extension of sensitive detector layer relevant for seeding in
  // negative direction in z in mm
  float zMin;
  // maximum distance in r from middle space point to bottom or top spacepoint
  // in mm
  float deltaRMax;
  // maximum forward direction expressed as cot(theta)
  float cotThetaMax;
  // maximum impact parameter in mm
  float impactMax;
  // number of consecutive phi bins used in the seeding
  int numberOfPhiBins = 1;
  // enable non equidistant binning in z
  std::vector < float > zBinEdges;

};


template <typename external_spacepoint_t>
using SpacePointGrid = detail::Grid<std::vector<std::unique_ptr<
                     const InternalSpacePoint<external_spacepoint_t>>>,
                 detail::Axis<detail::AxisType::Equidistant,
                              detail::AxisBoundaryType::Closed>,
                 detail::Axis<detail::AxisType::Variable,
                              detail::AxisBoundaryType::Bound>>;


class SpacePointGridCreator {
 public:
  template <typename external_spacepoint_t>
  static std::unique_ptr<SpacePointGrid<external_spacepoint_t>> createGrid(
      const Acts::SpacePointGridConfig& config);
};
}  // namespace Acts
#include "Acts/Seeding/SpacePointGrid.ipp"
