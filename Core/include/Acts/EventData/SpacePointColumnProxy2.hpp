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
  using IndexRange = SpacePointIndexRange2;
  using IndexSubset = SpacePointIndexSubset2;
  using Value = T;
  using Container = const_if_t<ReadOnly, SpacePointContainer2>;
  using Column = const_if_t<ReadOnly, std::vector<Value>>;

  /// Constructs a space point column proxy for the given container and column.
  /// @param container The container holding the space point.
  /// @param column The column of data to access.
  SpacePointColumnProxy(Container &container, Column &column)
      : m_container{&container}, m_column(&column) {}

  /// Copy construct a space point column proxy.
  /// @param other The space point column proxy to copy.
  SpacePointColumnProxy(const SpacePointColumnProxy &other) noexcept = default;

  /// Copy construct a mutable space point column proxy.
  /// @param other The mutable space point column proxy to copy.
  explicit SpacePointColumnProxy(
      const SpacePointColumnProxy<T, false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_column(&other.column()) {}

  /// Returns a const proxy of the space point column.
  /// @return A const proxy of the space point column.
  SpacePointColumnProxy<T, true> asConst() const noexcept
    requires(!ReadOnly)
  {
    return {*m_container, *m_column};
  }

  /// Gets the container holding the space point.
  /// @return A reference to the container holding the space point.
  SpacePointContainer2 &container() noexcept
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
  /// @param index The index of the space point to access.
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
  /// @param index The index of the space point to access.
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
  /// @param index The index of the space point to access.
  /// @return A mutable reference to the column entry at the given index.
  Value &operator[](Index index) noexcept
    requires(!ReadOnly)
  {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }
  /// Returns a const reference to the column entry at the given index.
  /// @param index The index of the space point to access.
  /// @return A const reference to the column entry at the given index.
  const Value &operator[](Index index) const noexcept {
    assert(index < column().size() && "Index out of bounds");
    return data()[index];
  }

  class Subset {
   public:
    class Iterator {
     public:
      using SubsetIterator = IndexSubset::iterator;

      using value_type = T;
      using difference_type = std::ptrdiff_t;

      using iterator_category = std::random_access_iterator_tag;
      using iterator_concept = std::random_access_iterator_tag;

      constexpr Iterator() noexcept = default;
      constexpr Iterator(Column &column, SubsetIterator iterator) noexcept
          : m_column(&column), m_iterator(iterator) {}

      constexpr value_type operator*() const noexcept {
        return (*m_column)[*m_iterator];
      }
      constexpr value_type operator[](difference_type n) const noexcept {
        return (*m_column)[m_iterator[n]];
      }

      constexpr Iterator &operator++() noexcept {
        ++m_iterator;
        return *this;
      }
      constexpr Iterator operator++(int) noexcept {
        auto tmp = *this;
        ++(*this);
        return tmp;
      }
      constexpr Iterator &operator--() noexcept {
        --m_iterator;
        return *this;
      }
      constexpr Iterator operator--(int) noexcept {
        auto tmp = *this;
        --(*this);
        return tmp;
      }

      constexpr Iterator &operator+=(difference_type n) noexcept {
        m_iterator += n;
        return *this;
      }
      constexpr Iterator &operator-=(difference_type n) noexcept {
        m_iterator -= n;
        return *this;
      }

     private:
      Column *m_column{};
      SubsetIterator m_iterator{};

      friend constexpr Iterator operator+(Iterator it,
                                          difference_type n) noexcept {
        return it += n;
      }

      friend constexpr Iterator operator+(difference_type n,
                                          Iterator it) noexcept {
        return it += n;
      }

      friend constexpr Iterator operator-(Iterator it,
                                          difference_type n) noexcept {
        return it -= n;
      }

      friend constexpr difference_type operator-(const Iterator &lhs,
                                                 const Iterator &rhs) noexcept {
        return lhs.m_iterator - rhs.m_iterator;
      }

      friend constexpr auto operator<=>(const Iterator &a,
                                        const Iterator &b) noexcept {
        return a.m_iterator <=> b.m_iterator;
      }
      friend constexpr bool operator==(const Iterator &a,
                                       const Iterator &b) noexcept {
        return a.m_iterator == b.m_iterator;
      }
    };
    using iterator = Iterator;

    constexpr Subset(Column &column, const IndexSubset &subset) noexcept
        : m_column(&column), m_subset(subset) {}

    constexpr Container &column() const noexcept { return *m_column; }
    constexpr const IndexSubset &subset() const noexcept { return m_subset; }

    constexpr std::size_t size() const noexcept { return m_subset.size(); }
    constexpr bool empty() const noexcept { return size() == 0; }

    constexpr auto front() const noexcept { return column()[m_subset.front()]; }
    constexpr auto back() const noexcept { return column()[m_subset.back()]; }

    constexpr iterator begin() const noexcept {
      return iterator(*m_column, m_subset.begin());
    }
    constexpr iterator end() const noexcept {
      return iterator(*m_column, m_subset.end());
    }

   private:
    Column *m_column{};
    IndexSubset m_subset{};
  };

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

template <typename T>
using ConstSpacePointColumnProxy = SpacePointColumnProxy<T, true>;
template <typename T>
using MutableSpacePointColumnProxy = SpacePointColumnProxy<T, false>;

}  // namespace Acts::Experimental
