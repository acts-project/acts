// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/ContainerSubset.hpp"

#include <cassert>
#include <span>
#include <vector>

namespace Acts {

class SpacePointContainer2;

/// Additional column of data that can be added to the spacepoint container.
/// The column is indexed by the spacepoint index.
template <typename T, bool read_only>
class SpacePointColumnProxy {
 public:
  /// Flag indicating whether this spacepoint column proxy is read-only
  constexpr static bool ReadOnly = read_only;
  /// Type alias for spacepoint index type
  using Index = SpacePointIndex2;
  /// Type alias for spacepoint index range type
  using IndexRange = SpacePointIndexRange2;
  /// Type alias for spacepoint index subset type
  using IndexSubset = SpacePointIndexSubset2;
  /// Type alias for column value type
  using Value = T;
  /// Type alias for container type (const if read-only)
  using Container = const_if_t<ReadOnly, SpacePointContainer2>;
  /// Type alias for column container type (const if read-only)
  using Column = const_if_t<ReadOnly, std::vector<Value>>;

  /// Constructs a spacepoint column proxy for the given container and column.
  /// @param container The container holding the spacepoint.
  /// @param column The column of data to access.
  SpacePointColumnProxy(Container &container, Column &column)
      : m_container{&container}, m_column(&column) {}

  /// Copy construct a spacepoint column proxy.
  /// @param other The spacepoint column proxy to copy.
  SpacePointColumnProxy(const SpacePointColumnProxy &other) noexcept = default;

  /// Copy construct a mutable spacepoint column proxy.
  /// @param other The mutable spacepoint column proxy to copy.
  explicit SpacePointColumnProxy(
      const SpacePointColumnProxy<T, false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_column(&other.column()) {}

  /// Returns a const proxy of the spacepoint column.
  /// @return A const proxy of the spacepoint column.
  SpacePointColumnProxy<T, true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, *m_column};
  }

  /// Gets the container holding the spacepoint.
  /// @return A reference to the container holding the spacepoint.
  SpacePointContainer2 &container() noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the spacepoint.
  /// @return A const reference to the container holding the spacepoint.
  const SpacePointContainer2 &container() const noexcept {
    return *m_container;
  }

  /// Returns a const reference to the column container.
  /// @return A const reference to the column container.
  const std::vector<Value> &column() const noexcept { return *m_column; }

  /// Returns a mutable span to the column data.
  /// @return A mutable span to the column data.
  std::span<Value> data() noexcept
    requires(!ReadOnly)
  {
    return std::span<Value>(column().data(), column().size());
  }
  /// Returns a const span to the column data.
  /// @return A const span to the column data.
  std::span<const Value> data() const noexcept {
    return std::span<const Value>(column().data(), column().size());
  }

  /// Returns a mutable reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the spacepoint to access.
  /// @return A mutable reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  Value &at(Index index)
    requires(!ReadOnly)
  {
    if (index >= column().size()) {
      throw std::out_of_range("Index out of range in SpacePointContainer2: " +
                              std::to_string(index) +
                              " >= " + std::to_string(size()));
    }
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the spacepoint to access.
  /// @return A const reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  const Value &at(Index index) const {
    if (index >= column().size()) {
      throw std::out_of_range("Index out of range in SpacePointContainer2: " +
                              std::to_string(index) +
                              " >= " + std::to_string(size()));
    }
    return data()[index];
  }

  /// Returns a mutable reference to the column entry at the given index.
  /// @param index The index of the spacepoint to access.
  /// @return A mutable reference to the column entry at the given index.
  Value &operator[](Index index) noexcept
    requires(!ReadOnly)
  {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// @param index The index of the spacepoint to access.
  /// @return A const reference to the column entry at the given index.
  const Value &operator[](Index index) const noexcept {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }

  /// Subset view over selected column entries.
  class Subset : public detail::ContainerSubset<Subset, Subset, Column, Value,
                                                Index, ReadOnly> {
   public:
    /// Base class type
    using Base =
        detail::ContainerSubset<Subset, Subset, Column, Value, Index, ReadOnly>;

    using Base::Base;
  };

  /// Creates a subset view of this spacepoint column based on provided
  /// indices.
  ///
  /// This method creates a subset proxy that provides access to only the space
  /// points at the indices specified in the IndexSubset. The subset maintains a
  /// reference to the original column data without copying, enabling efficient
  /// access to selected spacepoints for filtering, clustering, or other
  /// operations.
  ///
  /// @param subset The index subset specifying which spacepoints to include
  /// @return A subset proxy providing access to the selected spacepoints
  ///
  /// @note The returned subset shares data with the original column
  /// @note The subset remains valid only as long as the original column exists
  /// @note This operation does not copy data, providing efficient subset access
  Subset subset(const IndexSubset &subset) const noexcept {
    return Subset(*m_column, subset);
  }

 private:
  Container *m_container{};
  Column *m_column{};

  std::uint32_t size() const noexcept { return column().size(); }

  std::vector<Value> &column() noexcept
    requires(!ReadOnly)
  {
    return *m_column;
  }

  friend class SpacePointContainer2;
};

/// Const proxy to a spacepoint column for read-only access
template <typename T>
using ConstSpacePointColumnProxy = SpacePointColumnProxy<T, true>;
/// Mutable proxy to a spacepoint column allowing modification
template <typename T>
using MutableSpacePointColumnProxy = SpacePointColumnProxy<T, false>;

}  // namespace Acts
