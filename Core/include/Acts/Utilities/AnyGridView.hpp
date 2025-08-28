// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <type_traits>

namespace Acts {

/// @brief Type-safe view into a grid with mutable access
///
/// @tparam T Type of values stored in the grid
///
/// This class provides a type-safe interface to access grid objects through the
/// type-erased IGrid interface with mutable access. It ensures that the grid
/// being viewed contains values of the expected type T.
///
/// Example usage:
/// ```
/// Grid<double, Axis> grid(...);
/// AnyGridView<double> view(grid);
/// view.atLocalBins({1}) = 42.0;  // Modify the grid through the view
/// ```
template <typename T>
class AnyGridView {
 public:
  /// Type of pointer to grid, const or non-const depending on isConst
  using GridPointerType = IGrid*;
  /// Type for indices, imported from IGrid
  using AnyIndexType = IGrid::AnyIndexType;
  /// Type for points, imported from IGrid
  using AnyPointType = IGrid::AnyPointType;

  /// @brief Constructor from non-const IGrid reference
  /// @param grid The grid to view
  explicit AnyGridView(IGrid& grid) : m_grid(&grid) { checkType(); }

  /// @brief Constructor from non-const concrete Grid reference
  /// @tparam Axes Parameter pack of axis types defining the grid
  /// @param grid The concrete grid to view
  /// @note This constructor creates a mutable view of the grid
  template <typename... Axes>
  explicit AnyGridView(Grid<T, Axes...>& grid) : m_grid(&grid) {
    checkType();
  }

  /// Copy constructor
  /// @param other The AnyGridView to copy from
  AnyGridView(const AnyGridView& other) = default;
  /// Copy assignment operator
  /// @param other The AnyGridView to copy from
  /// @return Reference to this AnyGridView after copying
  AnyGridView& operator=(const AnyGridView& other) = default;

  /// Move constructor
  AnyGridView(AnyGridView&&) noexcept = default;
  /// Move assignment operator
  /// @return Reference to this AnyGridView after moving
  AnyGridView& operator=(AnyGridView&&) noexcept = default;

  /// @brief Access value at given local bin indices with mutable access
  /// @param indices The local bin indices
  /// @return Reference to the value at the specified bin
  /// @note This method provides mutable access to grid values
  /// @throws std::invalid_argument if indices size doesn't match grid dimensions
  /// @throws std::out_of_range if indices are out of bounds
  T& atLocalBins(const AnyIndexType& indices) {
    std::any any = m_grid->atLocalBinsAny(indices);
    return *std::any_cast<T*>(any);
  }

  /// @brief Access value at given local bin indices with const access
  /// @param indices The local bin indices
  /// @return Const reference to the value at the specified bin
  /// @throws std::invalid_argument if indices size doesn't match grid dimensions
  /// @throws std::out_of_range if indices are out of bounds
  const T& atLocalBins(const AnyIndexType& indices) const {
    std::any any = m_grid->atLocalBinsAny(indices);
    return *std::any_cast<T const*>(any);
  }

  /// @brief Get the number of dimensions of the grid
  /// @return The number of dimensions
  std::size_t dimensions() const { return m_grid->dimensions(); }

  /// @brief Get the center position of a bin for given indices
  /// @param indices The local bin indices
  /// @return The center position of the bin
  AnyPointType binCenter(const IGrid::AnyIndexType& indices) const {
    return m_grid->binCenterAny(indices);
  }

  /// @brief Get the lower left edge position of a bin for given indices
  /// @param indices The local bin indices
  /// @return The lower left edge position of the bin
  AnyPointType lowerLeftBinEdge(const IGrid::AnyIndexType& indices) const {
    return m_grid->lowerLeftBinEdgeAny(indices);
  }

  /// @brief Get the upper right edge position of a bin for given indices
  /// @param indices The local bin indices
  /// @return The upper right edge position of the bin
  AnyPointType upperRightBinEdge(const IGrid::AnyIndexType& indices) const {
    return m_grid->upperRightBinEdgeAny(indices);
  }

  /// @brief Get the number of bins along each axis
  /// @return Vector containing the number of bins for each axis
  AnyIndexType numLocalBins() const { return m_grid->numLocalBinsAny(); }

  /// @brief Check if the grid's value type matches the template parameter T
  /// @throws std::invalid_argument if there's a type mismatch
  void checkType() {
    if (m_grid->valueType() != typeid(T)) {
      throw std::invalid_argument("Type mismatch between grid and view type");
    }
  }

  /// Type-erased pointer to the underlying grid
  GridPointerType m_grid;
};

/// @brief Deduction guide for AnyGridView from Grid
/// @tparam T Type of values stored in the grid
/// @tparam Axes Parameter pack of axis types defining the grid
/// @param grid The grid to create a view for
template <typename T, typename... Axes>
AnyGridView(Grid<T, Axes...>& grid) -> AnyGridView<T>;

/// @brief Type-safe view into a grid with const access
///
/// @tparam T Type of values stored in the grid
///
/// This class provides a type-safe interface to access grid objects through the
/// type-erased IGrid interface with const access only. It ensures that the grid
/// being viewed contains values of the expected type T.
///
/// Example usage:
/// ```
/// const Grid<double, Axis> grid(...);
/// AnyGridConstView<double> view(grid);
/// double value = view.atLocalBins({1});  // Read-only access
/// ```
template <typename T>
class AnyGridConstView {
 public:
  /// Type of pointer to grid, const or non-const depending on isConst
  using GridPointerType = const IGrid*;
  /// Type for indices, imported from IGrid
  using AnyIndexType = IGrid::AnyIndexType;
  /// Type for points, imported from IGrid
  using AnyPointType = IGrid::AnyPointType;

  /// @brief Constructor from const IGrid reference
  /// @param grid The grid to view
  /// @note This constructor creates a const view of the grid
  explicit AnyGridConstView(const IGrid& grid) : m_grid(&grid) { checkType(); }

  /// @brief Constructor from const concrete Grid reference
  /// @tparam Axes Parameter pack of axis types defining the grid
  /// @param grid The concrete grid to view
  /// @note This constructor creates a const view of the grid
  template <typename... Axes>
  explicit AnyGridConstView(const Grid<T, Axes...>& grid) : m_grid(&grid) {
    checkType();
  }

  /// @brief Access value at given local bin indices with const access
  /// @param indices The local bin indices
  /// @return Const reference to the value at the specified bin
  /// @throws std::invalid_argument if indices size doesn't match grid dimensions
  /// @throws std::out_of_range if indices are out of bounds
  const T& atLocalBins(const AnyIndexType& indices) const {
    std::any any = m_grid->atLocalBinsAny(indices);
    return *std::any_cast<T const*>(any);
  }

  /// @brief Get the number of dimensions of the grid
  /// @return The number of dimensions
  std::size_t dimensions() const { return m_grid->dimensions(); }

  /// @brief Get the center position of a bin for given indices
  /// @param indices The local bin indices
  /// @return The center position of the bin
  AnyPointType binCenter(const IGrid::AnyIndexType& indices) const {
    return m_grid->binCenterAny(indices);
  }

  /// @brief Get the lower left edge position of a bin for given indices
  /// @param indices The local bin indices
  /// @return The lower left edge position of the bin
  AnyPointType lowerLeftBinEdge(const IGrid::AnyIndexType& indices) const {
    return m_grid->lowerLeftBinEdgeAny(indices);
  }

  /// @brief Get the upper right edge position of a bin for given indices
  /// @param indices The local bin indices
  /// @return The upper right edge position of the bin
  AnyPointType upperRightBinEdge(const IGrid::AnyIndexType& indices) const {
    return m_grid->upperRightBinEdgeAny(indices);
  }

  /// @brief Get the number of bins along each axis
  /// @return Vector containing the number of bins for each axis
  AnyIndexType numLocalBins() const { return m_grid->numLocalBinsAny(); }

  /// @brief Check if the grid's value type matches the template parameter T
  /// @throws std::invalid_argument if there's a type mismatch
  void checkType() {
    if (m_grid->valueType() != typeid(T)) {
      throw std::invalid_argument("Type mismatch between grid and view type");
    }
  }

  /// Type-erased pointer to the underlying const grid
  GridPointerType m_grid;
};

/// @brief Deduction guide for AnyGridConstView from const Grid
/// @tparam T Type of values stored in the grid
/// @tparam Axes Parameter pack of axis types defining the grid
/// @param grid The const grid to create a view for
template <typename T, typename... Axes>
AnyGridConstView(const Grid<T, Axes...>& grid) -> AnyGridConstView<T>;

}  // namespace Acts
