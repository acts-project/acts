// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <numbers>
#include <vector>

namespace Acts {

/// A cylindrical space point grid used for seeding in a cylindrical detector
/// geometry.
/// The grid is defined in cylindrical coordinates (phi, z, r) and allows for
/// efficient access to space points based on their azimuthal angle,
/// z-coordinate, and radial distance.
class CylindricalSpacePointGrid2 {
 public:
  /// Space point index type used in the grid.
  using SpacePointIndex = std::uint32_t;
  /// Type alias for bin container holding space point indices
  using BinType = std::vector<SpacePointIndex>;
  /// Type alias for phi axis with equidistant binning and closed boundaries
  using PhiAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Closed>;
  /// Type alias for z axis with variable binning and open boundaries
  using ZAxisType = Axis<AxisType::Variable, AxisBoundaryType::Open>;
  /// Type alias for r axis with variable binning and open boundaries
  using RAxisType = Axis<AxisType::Variable, AxisBoundaryType::Open>;
  /// Cylindrical space point grid type, which is a grid over `BinType` with
  /// axes defined by `PhiAxisType`, `ZAxisType`, and `RAxisType`.
  /// The grid is a 3D grid with the axes representing azimuthal angle (phi),
  /// z-coordinate, and radial distance (r).
  using GridType = Grid<BinType, PhiAxisType, ZAxisType, RAxisType>;
  /// Type alias for binned group over the cylindrical grid
  using BinnedGroupType = BinnedGroup<GridType>;

  /// Configuration parameters for the cylindrical space point grid.
  struct Config {
    /// minimum pT
    float minPt = 0 * UnitConstants::MeV;
    /// minimum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    /// WARNING: if rMin is smaller than impactMax, the bin size will be 2*pi,
    /// which will make seeding very slow!
    float rMin = 0 * UnitConstants::mm;
    /// maximum extension of sensitive detector layer relevant for seeding as
    /// distance from x=y=0 (i.e. in r)
    float rMax = 600 * UnitConstants::mm;
    /// minimum extension of sensitive detector layer relevant for seeding in
    /// negative direction in z
    float zMin = -2800 * UnitConstants::mm;
    /// maximum extension of sensitive detector layer relevant for seeding in
    /// positive direction in z
    float zMax = 2800 * UnitConstants::mm;
    /// maximum distance in r from middle space point to bottom or top
    /// spacepoint
    float deltaRMax = 270 * UnitConstants::mm;
    /// maximum forward direction expressed as cot(theta)
    float cotThetaMax = 10.01788;  // equivalent to eta = 3 (pseudorapidity)
    /// maximum impact parameter in mm
    float impactMax = 0 * UnitConstants::mm;
    /// minimum phi value for phiAxis construction
    float phiMin = -std::numbers::pi_v<float>;
    /// maximum phi value for phiAxis construction
    float phiMax = std::numbers::pi_v<float>;
    /// Multiplicator for the number of phi-bins. The minimum number of phi-bins
    /// depends on min_pt, magnetic field: 2*pi/(minPT particle phi-deflection).
    /// phiBinDeflectionCoverage is a multiplier for this number. If
    /// numPhiNeighbors (in the configuration of the BinFinders) is configured
    /// to return 1 neighbor on either side of the current phi-bin (and you want
    /// to cover the full phi-range of minPT), leave this at 1.
    int phiBinDeflectionCoverage = 1;
    /// maximum number of phi bins
    int maxPhiBins = 10000;
    /// enable non equidistant binning in z
    std::vector<float> zBinEdges{};
    /// enable non equidistant binning in r
    std::vector<float> rBinEdges{};

    /// magnetic field
    float bFieldInZ = 0 * UnitConstants::T;

    /// bin finder for bottom space points
    std::optional<GridBinFinder<3ul>> bottomBinFinder;
    /// bin finder for top space points
    std::optional<GridBinFinder<3ul>> topBinFinder;
    /// navigation structure for the grid
    std::array<std::vector<std::size_t>, 3ul> navigation;
  };

  /// Construct a cylindrical space point grid with the given configuration and
  /// an optional logger.
  /// @param config Configuration for the cylindrical grid
  /// @param logger Optional logger instance for debugging output
  explicit CylindricalSpacePointGrid2(
      const Config& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("CylindricalSpacePointGrid2", Logging::Level::INFO));

  /// Clear the grid and drop all state. The object will behave like a newly
  /// constructed one.
  void clear();

  /// Get the bin index for a space point given its azimuthal angle, radial
  /// distance, and z-coordinate.
  /// @param position The position of the space point in (phi, z, r) coordinates
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(const Vector3& position) const {
    if (!grid().isInside(position)) {
      return std::nullopt;
    }
    return grid().globalBinFromPosition(position);
  }
  /// Get the bin index for a space point given its azimuthal angle, radial
  /// distance, and z-coordinate.
  /// @param phi The azimuthal angle of the space point in radians
  /// @param z The z-coordinate of the space point
  /// @param r The radial distance of the space point from the origin
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(float phi, float z, float r) const {
    return binIndex(Vector3(phi, z, r));
  }

  /// Insert a space point into the grid.
  /// @param index The index of the space point to insert
  /// @param phi The azimuthal angle of the space point in radians
  /// @param z The z-coordinate of the space point
  /// @param r The radial distance of the space point from the origin
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(SpacePointIndex index, float phi, float z,
                                    float r);
  /// Insert a space point into the grid.
  /// @param sp The space point to insert
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(const ConstSpacePointProxy2& sp) {
    return insert(sp.index(), sp.phi(), sp.z(), sp.r());
  }

  /// Fill the grid with space points from the container.
  /// @param spacePoints The space point container to fill the grid with
  void extend(const SpacePointContainer2::ConstRange& spacePoints);

  /// Sort the bins in the grid by the space point radius, which is required by
  /// some algorithms that operate on the grid.
  /// @param spacePoints The space point container to sort the bins by radius
  void sortBinsByR(const SpacePointContainer2& spacePoints);

  /// Compute the range of radii in the grid. This requires the grid to be
  /// filled with space points and sorted by radius. The sorting can be done
  /// with the `sortBinsByR` method.
  /// @param spacePoints The space point container to compute the radius range
  /// @return The range of radii in the grid
  Range1D<float> computeRadiusRange(
      const SpacePointContainer2& spacePoints) const;

  /// Mutable bin access by index.
  /// @param index The index of the bin to access
  /// @return Mutable reference to the bin at the specified index
  BinType& at(std::size_t index) { return grid().at(index); }
  /// Const bin access by index.
  /// @param index The index of the bin to access
  /// @return Const reference to the bin at the specified index
  const BinType& at(std::size_t index) const { return grid().at(index); }

  /// Mutable grid access.
  /// @return Mutable reference to the grid
  GridType& grid() { return *m_grid; }
  /// Const grid access.
  /// @return Const reference to the grid
  const GridType& grid() const { return *m_grid; }

  /// Access to the binned group.
  /// @return Reference to the binned group
  const BinnedGroupType& binnedGroup() const { return *m_binnedGroup; }

  /// Get the number of space points in the grid.
  /// @return The number of space points in the grid
  std::size_t numberOfSpacePoints() const { return m_counter; }

  /// Get the number of bins in the grid.
  /// @return The number of bins in the grid
  std::size_t numberOfBins() const { return grid().size(); }

 private:
  Config m_cfg;

  std::unique_ptr<const Logger> m_logger;

  GridType* m_grid{};
  std::optional<BinnedGroupType> m_binnedGroup;

  std::size_t m_counter{};

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
