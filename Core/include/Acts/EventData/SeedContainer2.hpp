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
  std::size_t size() const noexcept { return m_sizes.size(); }
  /// Checks if the seed container is empty.
  /// @return True if the container is empty, false otherwise.
  bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of seeds.
  /// @param size The number of seeds to reserve space for.
  /// @param averageSpacePoints The average number of space points per seed.
  void reserve(std::size_t size, float averageSpacePoints = 3) noexcept;

  /// Clears the seed container, removing all seeds and space points.
  void clear() noexcept;

  /// Creates a new seed with the given space points.
  /// @param spacePoints The space points that make up the seed.
  /// @return A mutable proxy to the newly created seed.
  MutableProxyType createSeed(
      std::span<const SpacePointIndex2> spacePoints) noexcept;

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
  std::span<SpacePointIndex2> spacePointIndices(IndexType index) noexcept {
    assert(index < m_sizes.size() && "Index out of bounds");
    assert(index < m_spacePointOffsets.size() && "Index out of bounds");
    return std::span<SpacePointIndex2>(
        m_spacePoints.data() + m_spacePointOffsets[index],
        m_spacePoints.data() + m_spacePointOffsets[index] + m_sizes[index]);
  }

  /// Mutable access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the quality of the seed at the given index.
  float &quality(IndexType index) noexcept {
    assert(index < m_qualities.size() && "Index out of bounds");
    return m_qualities[index];
  }
  /// Mutable access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed at the
  float &vertexZ(IndexType index) noexcept {
    assert(index < m_vertexZs.size() && "Index out of bounds");
    return m_vertexZs[index];
  }

  /// Const access to the space point indices of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A span of space point indices associated with the seed at the given
  ///         index.
  std::span<const SpacePointIndex2> spacePointIndices(
      IndexType index) const noexcept {
    assert(index < m_sizes.size() && "Index out of bounds");
    assert(index < m_spacePointOffsets.size() && "Index out of bounds");
    return std::span<const SpacePointIndex2>(
        m_spacePoints.data() + m_spacePointOffsets[index],
        m_spacePoints.data() + m_spacePointOffsets[index] + m_sizes[index]);
  }

  /// Const access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the quality of the seed at the given index.
  float quality(IndexType index) const noexcept {
    assert(index < m_qualities.size() && "Index out of bounds");
    return m_qualities[index];
  }
  /// Const access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the vertex Z coordinate of the seed at the
  ///         given index.
  float vertexZ(IndexType index) const noexcept {
    assert(index < m_vertexZs.size() && "Index out of bounds");
    return m_vertexZs[index];
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

    value_type operator*() const { return value_type(*m_container, m_index); }

   private:
    ContainerType *m_container{};
    IndexType m_index{};

    friend bool operator==(const SeedIterator &a, const SeedIterator &b) {
      return a.m_index == b.m_index && a.m_container == b.m_container;
    }
  };
  using iterator = SeedIterator<false>;
  using const_iterator = SeedIterator<true>;

  iterator begin() noexcept { return iterator(*this, 0); }
  iterator end() noexcept { return iterator(*this, size()); }

  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

 private:
  std::vector<std::uint8_t> m_sizes;
  std::vector<std::size_t> m_spacePointOffsets;
  std::vector<float> m_qualities;
  std::vector<float> m_vertexZs;
  std::vector<SpacePointIndex2> m_spacePoints;
};

}  // namespace Acts::Experimental

#include "Acts/EventData/SeedContainer2.ipp"
