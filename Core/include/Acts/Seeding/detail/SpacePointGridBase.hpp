// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

namespace Acts::detail {

/// CRTP base class for the space point grids used in seeding. It owns the
/// grid and the binned group and implements the functionality that is
/// independent of the grid coordinate system. The derived class constructs
/// the axes and hands the grid over via `initializeGrid`, and maps a space
/// point onto the grid axes by implementing
/// `insert(const ConstSpacePointProxy&)`, which is used by `extend`.
/// @tparam derived_t The deriving space point grid class
/// @tparam grid_t The underlying grid type
template <typename derived_t, typename grid_t>
class SpacePointGridBase {
 public:
  /// Space point index type used in the grid.
  using SpacePointIndex = std::uint32_t;
  /// Type alias for bin container holding space point indices
  using BinType = std::vector<SpacePointIndex>;
  /// Type of the underlying grid
  using GridType = grid_t;
  /// Type alias for binned group over the grid
  using BinnedGroupType = BinnedGroup<GridType>;

  /// Clear the grid and drop all state. The object will behave like a newly
  /// constructed one.
  void clear() {
    for (std::size_t i = 0; i < grid().size(); ++i) {
      grid().at(i).clear();
    }
    m_counter = 0;
  }

  /// Get the bin index for a position on the grid axes.
  /// @param position The position of the space point on the grid axes
  /// @return The index of the bin in which the space point is located, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> binIndex(const Vector3& position) const {
    if (!grid().isInside(position)) {
      return std::nullopt;
    }
    return grid().globalBinFromPosition(position);
  }

  /// Insert a space point into the grid.
  /// @param index The index of the space point to insert
  /// @param position The position of the space point on the grid axes
  /// @return The index of the bin in which the space point was inserted, or
  ///         `std::nullopt` if the space point is outside the grid bounds.
  std::optional<std::size_t> insert(SpacePointIndex index,
                                    const Vector3& position) {
    const std::optional<std::size_t> gridIndex = binIndex(position);
    if (gridIndex.has_value()) {
      grid().at(*gridIndex).push_back(index);
      ++m_counter;
    }
    return gridIndex;
  }

  /// Fill the grid with space points from the container.
  /// @param spacePoints The space point container to fill the grid with
  void extend(const SpacePointContainer::ConstRange& spacePoints) {
    ACTS_VERBOSE("Inserting " << spacePoints.size()
                              << " space points to the grid");

    for (const ConstSpacePointProxy& sp : spacePoints) {
      derived().insert(sp);
    }
  }

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

 protected:
  /// Construct the base with a logger. The grid has to be handed over with
  /// `initializeGrid` afterwards.
  /// @param logger Logger instance for debugging output
  explicit SpacePointGridBase(std::unique_ptr<const Logger> logger)
      : m_logger(std::move(logger)) {}

  ~SpacePointGridBase() = default;

  /// Take ownership of the fully constructed grid and set up the binned
  /// group. Has to be called exactly once by the derived constructor.
  /// @param grid The grid to take ownership of
  /// @param bottomBinFinder Bin finder for bottom space points
  /// @param topBinFinder Bin finder for top space points
  /// @param navigation Navigation structure for the grid
  void initializeGrid(
      GridType&& grid, const GridBinFinder<GridType::DIM>& bottomBinFinder,
      const GridBinFinder<GridType::DIM>& topBinFinder,
      std::array<std::vector<std::size_t>, GridType::DIM> navigation) {
    m_binnedGroup.emplace(std::move(grid), bottomBinFinder, topBinFinder,
                          std::move(navigation));
    m_grid = &m_binnedGroup->grid();
  }

  /// Sort the space points in each bin by the given projection.
  /// @param spacePoints The space point container the stored indices refer to
  /// @param projection Callable mapping a space point proxy to its sort key
  template <typename projection_t>
  void sortBinsBy(const SpacePointContainer& spacePoints,
                  const projection_t& projection) {
    ACTS_VERBOSE("Sorting the grid");

    for (std::size_t i = 0; i < grid().size(); ++i) {
      BinType& bin = grid().at(i);
      std::ranges::sort(bin, {}, [&](SpacePointIndex spIndex) {
        return projection(spacePoints[spIndex]);
      });
    }

    ACTS_VERBOSE(
        "Number of space points inserted (within grid range): " << m_counter);
  }

  /// Compute the range of the projection values in the grid. This requires
  /// the bins to be sorted by the same projection with `sortBinsBy`.
  /// @param spacePoints The space point container the stored indices refer to
  /// @param projection Callable mapping a space point proxy to its sort key
  /// @return The range of projection values in the grid
  template <typename projection_t>
  Range1D<float> computeRange(const SpacePointContainer& spacePoints,
                              const projection_t& projection) const {
    float minRange = std::numeric_limits<float>::max();
    float maxRange = std::numeric_limits<float>::lowest();
    for (const BinType& bin : grid()) {
      if (bin.empty()) {
        continue;
      }
      auto first = spacePoints[bin.front()];
      auto last = spacePoints[bin.back()];
      minRange = std::min(projection(first), minRange);
      maxRange = std::max(projection(last), maxRange);
    }
    return {minRange, maxRange};
  }

  /// Access to the logger.
  /// @return Reference to the logger
  const Logger& logger() const { return *m_logger; }

 private:
  derived_t& derived() { return static_cast<derived_t&>(*this); }

  std::unique_ptr<const Logger> m_logger;

  GridType* m_grid{};
  std::optional<BinnedGroupType> m_binnedGroup;

  std::size_t m_counter{};
};

}  // namespace Acts::detail
