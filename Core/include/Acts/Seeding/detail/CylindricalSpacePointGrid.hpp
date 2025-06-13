// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <numbers>
#include <vector>

namespace Acts {

/// Concept to check the provided external space point type
/// can be used to fill the space point grid
template <typename external_spacepoint_t>
concept CylindricalGridElement = requires(external_spacepoint_t sp) {
  { sp.phi() } -> std::same_as<float>;
  { sp.z() } -> std::same_as<float>;
  { sp.radius() } -> std::same_as<float>;
};

/// Cylindrical Space Point bin is a 2D grid with (phi, z) bins
/// It stores a vector of internal space points to external space points
template <CylindricalGridElement external_spacepoint_t>
using CylindricalSpacePointGrid =
    Grid<std::vector<const external_spacepoint_t*>,
         Axis<AxisType::Equidistant, AxisBoundaryType::Closed>,
         Axis<AxisType::Variable, AxisBoundaryType::Open>,
         Axis<AxisType::Variable, AxisBoundaryType::Open>>;

/// Cylindrical Binned Group
template <typename external_spacepoint_t>
using CylindricalBinnedGroup =
    BinnedGroup<CylindricalSpacePointGrid<external_spacepoint_t>>;

template <typename external_spacepoint_t>
using CylindricalBinnedGroupIterator =
    BinnedGroupIterator<CylindricalSpacePointGrid<external_spacepoint_t>>;

struct CylindricalSpacePointGridConfig {
  // minimum pT to be found by seedFinder
  float minPt = 0 * UnitConstants::MeV;
  // maximum extension of sensitive detector layer relevant for seeding as
  // distance from x=y=0 (i.e. in r)
  float rMax = 320 * UnitConstants::mm;
  // maximum extension of sensitive detector layer relevant for seeding as
  // distance from x=y=0 (i.e. in r)
  float rMin = 0 * UnitConstants::mm;
  // maximum extension of sensitive detector layer relevant for seeding in
  // positive direction in z
  float zMax = 0 * UnitConstants::mm;
  // maximum extension of sensitive detector layer relevant for seeding in
  // negative direction in z
  float zMin = 0 * UnitConstants::mm;
  // maximum distance in r from middle space point to bottom or top spacepoint
  float deltaRMax = 0 * UnitConstants::mm;
  // maximum forward direction expressed as cot(theta)
  float cotThetaMax = 0;
  // maximum impact parameter in mm
  float impactMax = 0 * UnitConstants::mm;
  // minimum phi value for phiAxis construction
  float phiMin = -std::numbers::pi_v<float>;
  // maximum phi value for phiAxis construction
  float phiMax = std::numbers::pi_v<float>;
  // Multiplicator for the number of phi-bins. The minimum number of phi-bins
  // depends on min_pt, magnetic field: 2*pi/(minPT particle phi-deflection).
  // phiBinDeflectionCoverage is a multiplier for this number. If
  // numPhiNeighbors (in the configuration of the BinFinders) is configured to
  // return 1 neighbor on either side of the current phi-bin (and you want to
  // cover the full phi-range of minPT), leave this at 1.
  int phiBinDeflectionCoverage = 1;
  // maximum number of phi bins
  int maxPhiBins = 10000;
  // enable non equidistant binning in z
  std::vector<float> zBinEdges{};
  std::vector<float> rBinEdges{};

  bool isInInternalUnits = true;
  //[[deprecated("CylindricalSpacePointGridConfig uses internal units")]]
  CylindricalSpacePointGridConfig toInternalUnits() const { return *this; }

  void checkConfig() const {
    if (phiMin < -std::numbers::pi_v<float> ||
        phiMax > std::numbers::pi_v<float>) {
      throw std::runtime_error(
          "CylindricalSpacePointGridConfig: phiMin (" + std::to_string(phiMin) +
          ") and/or phiMax (" + std::to_string(phiMax) +
          ") are outside the allowed phi range, defined as "
          "[-std::numbers::pi_v<float>, std::numbers::pi_v<float>]");
    }
    if (phiMin > phiMax) {
      throw std::runtime_error(
          "CylindricalSpacePointGridConfig: phiMin is bigger then phiMax");
    }
    if (rMin > rMax) {
      throw std::runtime_error(
          "CylindricalSpacePointGridConfig: rMin is bigger then rMax");
    }
    if (zMin > zMax) {
      throw std::runtime_error(
          "CylindricalSpacePointGridConfig: zMin is bigger than zMax");
    }
  }
};

struct CylindricalSpacePointGridOptions {
  // magnetic field
  float bFieldInZ = 0 * UnitConstants::T;

  bool isInInternalUnits = true;
  //[[deprecated("CylindricalSpacePointGridOptions uses internal units")]]
  CylindricalSpacePointGridOptions toInternalUnits() const { return *this; }
};

/// Instructions on how to create and fill this grid specialization
class CylindricalSpacePointGridCreator {
 public:
  template <typename external_spacepoint_t>
  static CylindricalSpacePointGrid<external_spacepoint_t> createGrid(
      const CylindricalSpacePointGridConfig& _config,
      const CylindricalSpacePointGridOptions& _options,
      const Logger& logger = getDummyLogger());

  template <typename external_spacepoint_t,
            typename external_spacepoint_iterator_t>
  static void fillGrid(const SeedFinderConfig<external_spacepoint_t>& config,
                       const SeedFinderOptions& options,
                       CylindricalSpacePointGrid<external_spacepoint_t>& grid,
                       external_spacepoint_iterator_t spBegin,
                       external_spacepoint_iterator_t spEnd,
                       const Logger& logger = getDummyLogger());

  template <typename external_spacepoint_t, typename external_collection_t>
    requires std::ranges::range<external_collection_t> &&
             std::same_as<typename external_collection_t::value_type,
                          external_spacepoint_t>
  static void fillGrid(const SeedFinderConfig<external_spacepoint_t>& config,
                       const SeedFinderOptions& options,
                       CylindricalSpacePointGrid<external_spacepoint_t>& grid,
                       const external_collection_t& collection,
                       const Logger& logger = getDummyLogger());
};

}  // namespace Acts

#include "Acts/Seeding/detail/CylindricalSpacePointGrid.ipp"
