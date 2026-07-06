// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <numbers>
#include <vector>

namespace Acts {

/// A cartesian space point grid used for seeding in a cartesian detector
/// geometry.
/// The grid is defined in cartesian coordinates (x,y,z).
class CartesianSpacePointGrid {
 public:
  /// Space point index type used in the grid.
  using SpacePointIndex = std::uint32_t;
  /// Type alias for bin container holding space point indices
  using BinType = std::vector<SpacePointIndex>;
  /// Type alias for x axis with equidistant binning and open boundaries
  using XAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  /// Type alias for y axis with equidistant binning and open boundaries
  using YAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  /// Type alias for z axis with equidistant binning and open boundaries
  using ZAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  /// cartesian space point grid type, which is a grid over `BinType` with
  /// axes defined by `XAxisType`, `YAxisType`, and `ZAxisType`.
  /// The grid is a 3D grid with the axes representing azimuthal angle (phi),
  /// z-coordinate, and radial distance (r).
  using GridType = Grid<BinType, XAxisType, YAxisType, ZAxisType>;
  /// Type alias for binned group over the cartesian grid
  using BinnedGroupType = BinnedGroup<GridType>;

  /// Configuration parameters for the cartesian space point grid.
  struct Config {
    /// minimum pT

    /// number of bins in the x-coordinate
    int nXbins = 10;
    /// minimum x-coordinate of space points used in the seeding
    float xMin = -600 * UnitConstants::mm;
    /// maximum x-coordinate of space points used in the seeding
    float xMax = 600 * UnitConstants::mm;
    /// number of bins in the y-coordinate
    int nYbins = 6;
    /// minimum y-coordinate of space points used in the seeding
    float yMin = -600 * UnitConstants::mm;
    /// maximum y-coordinate of space points used in the seeding
    float yMax = 600 * UnitConstants::mm;
    /// number of bins in the z-coordinate
    int nZbins = 10;
    /// minimum z coordinate ofspace points used in the seeding
    float zMin = -2700 * UnitConstants::mm;
    /// maximum z coordinate ofspace points used in the seeding
    float zMax = 2700 * UnitConstants::mm;

    /// bin finder for bottom space points
    std::optional<GridBinFinder<3ul>> bottomBinFinder;
    /// bin finder for top space points
    std::optional<GridBinFinder<3ul>> topBinFinder;
    /// navigation structure for the grid
    std::array<std::vector<std::size_t>, 3ul> navigation;

    /// coordinate to sort the space points by
    Acts::CoordinateIndices sortingCoord = Acts::CoordinateIndices::eY;

    /// direction to sort the space points in
    Acts::Direction sortingDirection = Acts::Direction::Negative();
  };

  /// Construct a cartesian space point grid with the given configuration and
  /// an optional logger.
  /// @param config Configuration for the cartesian grid
  /// @param logger Optional logger instance for debugging output
  explicit CartesianSpacePointGrid(
      const Config& config,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("CartesianSpacePointGrid", Logging::Level::INFO));

  /// Clear the grid and drop all state. The object will behave like a newly
  /// constructed one.
  void clear();

  /// Get the bin index for a space point given its x-, y-, and z-coordinate.
  /// @param position The position of the space point in (x,y,z) coordinates
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
  /// @param x The x-coordinate of the space point
  /// @param y The y-coordinate of the space point
  /// @param z The z-coordinate of the space point
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(float x, float y, float z) const {
    return binIndex(Vector3(x, y, z));
  }

  /// Insert a space point into the grid.
  /// @param index The index of the space point to insert
  /// @param x The x-coordinate the space point
  /// @param y The y-coordinate the space point
  /// @param z The z-coordinate the space point
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(SpacePointIndex index, float x, float y,
                                    float z);
  /// Insert a space point into the grid.
  /// @param sp The space point to insert
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(const ConstSpacePointProxy2& sp) {
    return insert(sp.index(), sp.x(), sp.y(), sp.z());
  }

  /// Fill the grid with space points from the container.
  /// @param spacePoints The space point container to fill the grid with
  void extend(const SpacePointContainer2::ConstRange& spacePoints);

  /// Sort the bins in the grid by the space point sorting-coordinate, which is
  /// required by some algorithms that operate on the grid.
  /// @param spacePoints The space point container to sort the bins by the configured coordinate
  void sortBinsByCoord(const SpacePointContainer2& spacePoints);

  /// Compute the range of the sorting coordinates in the grid. This requires
  /// the grid to be filled with space points and sorted by coordinate. The
  /// sorting can be done with the `sortBinsByY` method.
  /// @param spacePoints The space point container to compute the y-range
  /// @return The range of y-values in the grid
  Range1D<float> computeCoordRange(
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

  std::function<float(const ConstSpacePointProxy2&)> m_sortCoordGetter;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
