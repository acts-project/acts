// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/Types.hpp"

#include <cassert>
#include <limits>
#include <span>
#include <vector>

namespace Acts::Experimental {

template <bool read_only>
class SeedProxy2;

using MutableSeedProxy2 = SeedProxy2<false>;
using ConstSeedProxy2 = SeedProxy2<true>;

/// A container of seeds. Individual seeds are modeled as a sequence of N space
/// points which are addressed via an index into the space point container.
/// Individual seeds are addressed via index. A proxy object simplifies the
/// handling.
class SeedContainer2 {
 public:
  using IndexType = SeedIndex2;
  using MutableProxyType = MutableSeedProxy2;
  using ConstProxyType = ConstSeedProxy2;

  /// Returns the size of the seed container, i.e., the number of seeds
  /// contained in it.
  /// @return The number of seeds in the container.
  std::size_t size() const { return m_entries.size(); }
  /// Checks if the seed container is empty.
  /// @return True if the container is empty, false otherwise.
  bool empty() const { return size() == 0; }

  /// Reserves space for the given number of seeds.
  /// @param size The number of seeds to reserve space for.
  /// @param averageSpacePoints The average number of space points per seed.
  void reserve(std::size_t size, float averageSpacePoints = 3);

  /// Clears the seed container, removing all seeds and space points.
  void clear();

  /// Creates a new seed with the given space points.
  /// @param spacePoints The space points that make up the seed.
  /// @return A mutable proxy to the newly created seed.
  MutableProxyType createSeed(std::span<const SpacePointIndex2> spacePoints);

  /// Returns a mutable proxy to the seed at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the seed to access.
  /// @return A mutable proxy to the seed at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxyType at(IndexType index);
  /// Returns a const proxy to the seed at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the seed to access.
  /// @return A const proxy to the seed at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxyType at(IndexType index) const;

  /// Mutable access to the space point indices of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A span of space point indices associated with the seed at the given
  ///         index.
  std::span<SpacePointIndex2> spacePointIndices(IndexType index) {
    assert(index < m_entries.size() && "Index out of bounds");
    return std::span<SpacePointIndex2>(
        m_spacePoints.data() + m_entries[index].spacePointOffset,
        m_spacePoints.data() + m_entries[index].spacePointOffset +
            m_entries[index].seedSize);
  }
  /// Mutable access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the quality of the seed at the given index.
  float &quality(IndexType index) {
    assert(index < m_entries.size() && "Index out of bounds");
    return m_entries[index].quality;
  }
  /// Mutable access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed at the
  float &vertexZ(IndexType index) {
    assert(index < m_entries.size() && "Index out of bounds");
    return m_entries[index].vertexZ;
  }

  /// Const access to the space point indices of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A span of space point indices associated with the seed at the given
  ///         index.
  std::span<const SpacePointIndex2> spacePointIndices(IndexType index) const {
    assert(index < m_entries.size() && "Index out of bounds");
    return std::span<const SpacePointIndex2>(
        m_spacePoints.data() + m_entries[index].spacePointOffset,
        m_spacePoints.data() + m_entries[index].spacePointOffset +
            m_entries[index].seedSize);
  }
  /// Const access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the quality of the seed at the given index.
  float quality(IndexType index) const {
    assert(index < m_entries.size() && "Index out of bounds");
    return m_entries[index].quality;
  }
  /// Const access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the vertex Z coordinate of the seed at the
  ///         given index.
  float vertexZ(IndexType index) const {
    assert(index < m_entries.size() && "Index out of bounds");
    return m_entries[index].vertexZ;
  }

  template <bool read_only>
  class SeedIterator {
   public:
    static constexpr bool ReadOnly = read_only;

    using ContainerType = const_if_t<ReadOnly, SeedContainer2>;

    using iterator_category = std::forward_iterator_tag;
    using value_type = SeedProxy2<ReadOnly>;
    using difference_type = std::ptrdiff_t;

    SeedIterator() = default;
    SeedIterator(ContainerType &container, IndexType index)
        : m_container(&container), m_index(index) {}

    SeedIterator &operator++() {
      ++m_index;
      return *this;
    }
    SeedIterator operator++(int) {
      SeedIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    bool operator==(const SeedIterator &other) const {
      return m_index == other.m_index && m_container == other.m_container;
    }
    bool operator!=(const SeedIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};
  };
  using iterator = SeedIterator<false>;
  using const_iterator = SeedIterator<true>;

  iterator begin() { return iterator(*this, 0); }
  iterator end() { return iterator(*this, size()); }

  const_iterator begin() const { return const_iterator(*this, 0); }
  const_iterator end() const { return const_iterator(*this, size()); }

 private:
  struct Entry {
    std::size_t seedSize{};
    std::size_t spacePointOffset{};
    float quality{-std::numeric_limits<float>::infinity()};
    float vertexZ{};

    Entry() = default;
    Entry(std::size_t size, std::size_t offset, float q, float vz)
        : seedSize(size), spacePointOffset(offset), quality(q), vertexZ(vz) {}
  };

  std::vector<Entry> m_entries{};
  std::vector<SpacePointIndex2> m_spacePoints{};
};

}  // namespace Acts::Experimental

#include "Acts/EventData/SeedContainer2.ipp"
