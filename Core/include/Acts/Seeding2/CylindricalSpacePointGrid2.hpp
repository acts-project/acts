// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData2/SpacePointContainer2.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <numbers>
#include <vector>

namespace Acts {

class CylindricalSpacePointGrid2 {
 public:
  /// Cylindrical Space Point bin is a 2D grid with (phi, z) bins
  /// It stores a vector of internal space points to external space points
  using GridType = Grid<std::vector<SpacePointIndex2>,
                        Axis<AxisType::Equidistant, AxisBoundaryType::Closed>,
                        Axis<AxisType::Variable, AxisBoundaryType::Open>,
                        Axis<AxisType::Variable, AxisBoundaryType::Open>>;

  /// Cylindrical Binned Group
  using BinnedGroup = BinnedGroup<GridType>;

  using BinnedGroupIterator = BinnedGroupIterator<GridType>;

  struct DerivedConfig;

  struct Config {
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

    // magnetic field
    float bFieldInZ = 0. * UnitConstants::T;

    DerivedConfig derive() const;
  };

  struct DerivedConfig : public Config {};

  explicit CylindricalSpacePointGrid2(
      const DerivedConfig& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("CylindricalSpacePointGrid2", Logging::Level::INFO));

  void clear();

  void insert(const ConstSpacePointProxy2& sp) {
    // fill rbins into grid
    Vector3 position(sp.phi(), sp.z(), sp.radius());
    if (!grid().isInside(position)) {
      return;
    }

    std::size_t globIndex = grid().globalBinFromPosition(position);
    auto& rbin = grid().at(globIndex);
    rbin.push_back(sp.index());
    ++m_counter;

    // keep track of the bins we modify so that we can later sort the SPs in
    // those bins only
    if (rbin.size() > 1 && !m_usedBinIndex[globIndex]) {
      m_usedBinIndex[globIndex] = true;
      m_rBinsIndex.push_back(globIndex);
    }
  }

  void extend(const SpacePointContainer2::ConstRange& spacePoints);

  void fill(const SpacePointContainer2& spacePoints);

  void sort(const SpacePointContainer2& spacePoints);

  Range1D<float> computeRadiusRange(
      const SpacePointContainer2& spacePoints) const;

  GridType& grid() { return *m_grid; }
  const GridType& grid() const { return *m_grid; }

  BinnedGroup binnedGround(
      const Acts::GridBinFinder<3ul>& bottomBinFinder,
      const Acts::GridBinFinder<3ul>& topBinFinder,
      const std::array<std::vector<std::size_t>, 3ul>& navigation) && {
    return BinnedGroup(std::move(*m_grid), bottomBinFinder, topBinFinder,
                       navigation);
  }

 private:
  DerivedConfig m_cfg;

  std::unique_ptr<const Logger> m_logger;

  std::optional<GridType> m_grid;

  std::size_t m_counter{};
  std::vector<bool> m_usedBinIndex;
  std::vector<std::size_t> m_rBinsIndex;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
