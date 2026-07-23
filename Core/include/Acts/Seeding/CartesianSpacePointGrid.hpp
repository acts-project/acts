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
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/detail/SpacePointGridBase.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <functional>
#include <vector>

namespace Acts {

/// A cartesian space point grid used for seeding in a cartesian detector
/// geometry.
/// The grid is defined in cartesian coordinates (x,y,z).
class CartesianSpacePointGrid
    : public detail::SpacePointGridBase<
          CartesianSpacePointGrid,
          Grid<std::vector<std::uint32_t>,
               Axis<AxisType::Equidistant, AxisBoundaryType::Open>,
               Axis<AxisType::Equidistant, AxisBoundaryType::Open>,
               Axis<AxisType::Equidistant, AxisBoundaryType::Open>>> {
 public:
  /// Type alias for x axis with equidistant binning and open boundaries
  using XAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  /// Type alias for y axis with equidistant binning and open boundaries
  using YAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  /// Type alias for z axis with equidistant binning and open boundaries
  using ZAxisType = Axis<AxisType::Equidistant, AxisBoundaryType::Open>;
  /// Type alias for the CRTP base class
  using GridBase =
      detail::SpacePointGridBase<CartesianSpacePointGrid, GridType>;

  /// Configuration parameters for the cartesian space point grid.
  struct Config {
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

  using GridBase::binIndex;
  /// Get the bin index for a space point given its x-, y-, and z-coordinate.
  /// @param x The x-coordinate of the space point
  /// @param y The y-coordinate of the space point
  /// @param z The z-coordinate of the space point
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(float x, float y, float z) const {
    return binIndex(Vector3(x, y, z));
  }

  using GridBase::insert;
  /// Insert a space point into the grid.
  /// @param index The index of the space point to insert
  /// @param x The x-coordinate the space point
  /// @param y The y-coordinate the space point
  /// @param z The z-coordinate the space point
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(SpacePointIndex index, float x, float y,
                                    float z) {
    return insert(index, Vector3(x, y, z));
  }
  /// Insert a space point into the grid.
  /// @param sp The space point to insert
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(const ConstSpacePointProxy& sp) {
    return insert(sp.index(), sp.x(), sp.y(), sp.z());
  }

  /// Sort the bins in the grid by the space point sorting-coordinate, which is
  /// required by some algorithms that operate on the grid.
  /// @param spacePoints The space point container to sort the bins by the configured coordinate
  void sortBinsByCoord(const SpacePointContainer& spacePoints);

  /// Compute the range of the sorting coordinates in the grid. This requires
  /// the grid to be filled with space points and sorted by coordinate. The
  /// sorting can be done with the `sortBinsByCoord` method.
  /// @param spacePoints The space point container to compute the coordinate range
  /// @return The range of sorting-coordinate values in the grid
  Range1D<float> computeCoordRange(
      const SpacePointContainer& spacePoints) const;

 private:
  Config m_cfg;

  std::function<float(const ConstSpacePointProxy&)> m_sortCoordGetter;
};

}  // namespace Acts
