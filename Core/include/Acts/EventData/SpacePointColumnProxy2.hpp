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

#include <cassert>
#include <span>
#include <vector>

namespace Acts::Experimental {

class SpacePointContainer2;

/// Additional column of data that can be added to the space point container.
/// The column is indexed by the space point index.
template <typename T, bool read_only>
class SpacePointColumnProxy {
 public:
  constexpr static bool ReadOnly = read_only;
  using Index = SpacePointIndex2;
  using Value = const_if_t<ReadOnly, T>;
  using Container = const_if_t<ReadOnly, SpacePointContainer2>;
  using Column = const_if_t<ReadOnly, std::vector<T>>;
  using Span = std::span<Value>;
  using Reference = Value &;

  SpacePointColumnProxy(Container &container, Column &column)
      : m_container{&container}, m_column(&column) {}

  /// Gets the container holding the space point.
  /// @return A reference to the container holding the space point.
  SpacePointContainer2 &container() const noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the space point.
  /// @return A const reference to the container holding the space point.
  const SpacePointContainer2 &container() const noexcept {
    return *m_container;
  }

  /// Returns a const reference to the column container.
  /// @return A const reference to the column container.
  const Column &column() const noexcept { return *m_column; }

  /// Returns a span to the column data.
  /// @return A span to the column data.
  Span data() const noexcept { return Span(column().data(), column().size()); }

  /// Returns a reference to the column entry at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the space point to access.
  /// @return A reference to the column entry at the given index.
  /// @throws std::out_of_range if the index is out of range.
  Reference at(Index index) const {
    if (index >= column().size()) {
      throw std::out_of_range("Index out of range in SpacePointContainer2: " +
                              std::to_string(index) +
                              " >= " + std::to_string(size()));
    }
    return data()[index];
  }

  /// Returns a reference to the column entry at the given index.
  /// @param index The index of the space point to access.
  /// @return A reference to the column entry at the given index.
  Reference operator[](Index index) const noexcept {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }

 private:
  Container *m_container{};
  Column *m_column{};

  std::uint32_t size() const noexcept { return column().size(); }

  Column &column() const noexcept
    requires(!ReadOnly)
  {
    return *m_column;
  }

  friend class SpacePointContainer2;
};

template <typename T>
using ConstSpacePointColumnProxy = SpacePointColumnProxy<T, true>;
template <typename T>
using MutableSpacePointColumnProxy = SpacePointColumnProxy<T, false>;

}  // namespace Acts::Experimental
