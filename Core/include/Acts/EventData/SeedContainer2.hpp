// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"

#include <limits>
#include <span>

namespace Acts::Experimental {

using SeedIndex2 = std::size_t;

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
  void reserve(std::size_t size, float averageSpacePoints = 3) {
    m_entries.reserve(size);
    m_spacePoints.reserve(static_cast<std::size_t>(size * averageSpacePoints));
  }
  /// Clears the seed container, removing all seeds and space points.
  void clear() {
    m_entries.clear();
    m_spacePoints.clear();
  }

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
    return std::span<SpacePointIndex2>(
        m_spacePoints.data() + m_entries[index].spacePointOffset,
        m_spacePoints.data() + m_entries[index].spacePointOffset +
            m_entries[index].seedSize);
  }
  /// Mutable access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the quality of the seed at the given index.
  float &quality(IndexType index) { return m_entries[index].quality; }
  /// Mutable access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed at the
  float &vertexZ(IndexType index) { return m_entries[index].vertexZ; }

  /// Const access to the space point indices of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A span of space point indices associated with the seed at the given
  ///         index.
  std::span<const SpacePointIndex2> spacePointIndices(IndexType index) const {
    return std::span<const SpacePointIndex2>(
        m_spacePoints.data() + m_entries[index].spacePointOffset,
        m_spacePoints.data() + m_entries[index].spacePointOffset +
            m_entries[index].seedSize);
  }
  /// Const access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the quality of the seed at the given index.
  float quality(IndexType index) const { return m_entries[index].quality; }
  /// Const access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the vertex Z coordinate of the seed at the
  ///         given index.
  float vertexZ(IndexType index) const { return m_entries[index].vertexZ; }

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

/// A proxy class for accessing individual seeds.
template <bool read_only>
class SeedProxy2 {
 public:
  /// Indicates whether this seed proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  using IndexType = SeedIndex2;

  using ContainerType = const_if_t<ReadOnly, SeedContainer2>;

  /// Constructs a seed proxy for the given container and index.
  /// @param container The container holding the seed.
  /// @param index The index of the seed in the container.
  SeedProxy2(ContainerType &container, IndexType index)
      : m_container{&container}, m_index{index} {}

  /// Copy construct a seed proxy.
  /// @param other The seed proxy to copy.
  SeedProxy2(const SeedProxy2 &other) = default;

  /// Copy construct a mutable seed proxy.
  /// @param other The mutable seed proxy to copy.
  SeedProxy2(const SeedProxy2<false> &other)
    requires(ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  /// Gets the container holding the seed.
  /// @return A reference to the container holding the seed.
  SeedContainer2 &container()
    requires(ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the seed.
  /// @return A const reference to the container holding the seed.
  const SeedContainer2 &container() const { return *m_container; }
  /// Gets the index of the seed in the container.
  /// @return The index of the seed in the container.
  IndexType index() const { return m_index; }

  /// Returns the size of the seed, i.e., the number of space points
  /// associated with it.
  /// @return The number of space points in the seed.
  [[nodiscard]] std::size_t size() const {
    return m_container->spacePointIndices(m_index).size();
  }
  /// Checks if the seed is empty, i.e., has no space points associated with it.
  /// @return True if the seed is empty, false otherwise.
  [[nodiscard]]
  bool empty() const {
    return size() == 0;
  }

  /// Mutable access to the space point indices of the seed.
  /// @return A mutable span of space point indices associated with the seed.
  std::span<SpacePointIndex2> spacePointIndices()
    requires(!ReadOnly)
  {
    return m_container->spacePointIndices(m_index);
  }
  /// Mutable access to the quality of the seed.
  /// @return A mutable reference to the quality of the seed.
  float &quality()
    requires(!ReadOnly)
  {
    return m_container->quality(m_index);
  }
  /// Mutable access to the vertex Z coordinate of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed.
  float &vertexZ()
    requires(!ReadOnly)
  {
    return m_container->vertexZ(m_index);
  }

  /// Const access to the space point indices of the seed.
  /// @return A span of space point indices associated with the seed.
  ///         This span is read-only and cannot be modified.
  std::span<const SpacePointIndex2> spacePointIndices() const {
    return m_container->spacePointIndices(m_index);
  }
  /// Const access to the quality of the seed.
  /// @return The quality of the seed.
  float quality() const { return m_container->quality(m_index); }
  /// Const access to the vertex Z coordinate of the seed.
  /// @return The vertex Z coordinate of the seed.
  float vertexZ() const { return m_container->vertexZ(m_index); }

  class SpacePointIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = ConstSpacePointProxy2;
    using difference_type = std::ptrdiff_t;

    SpacePointIterator() = default;
    SpacePointIterator(const SpacePointContainer2 &spacePointContainer,
                       const SpacePointIndex2 *indexPointer)
        : m_spacePointContainer{&spacePointContainer},
          m_indexPointer{indexPointer} {}

    SpacePointIterator &operator++() {
      ++m_indexPointer;
      return *this;
    }

    SpacePointIterator operator++(int) {
      SpacePointIterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const SpacePointIterator &other) const {
      return m_indexPointer == other.m_indexPointer;
    }
    bool operator!=(const SpacePointIterator &other) const {
      return !(*this == other);
    }

    value_type operator*() const {
      return m_spacePointContainer->at(*m_indexPointer);
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    const SpacePointIndex2 *m_indexPointer{nullptr};
  };

  class SpacePointRange {
   public:
    SpacePointRange(const SpacePointContainer2 &spacePointContainer,
                    std::span<const SpacePointIndex2> spacePointIndices)
        : m_spacePointContainer{&spacePointContainer},
          m_spacePointIndices{spacePointIndices} {}

    std::size_t size() const { return m_spacePointIndices.size(); }
    bool empty() const { return size() == 0; }

    ConstSpacePointProxy2 operator[](std::size_t index) const {
      return m_spacePointContainer->at(m_spacePointIndices[index]);
    }

    SpacePointIterator begin() const {
      return SpacePointIterator(*m_spacePointContainer,
                                m_spacePointIndices.data());
    }
    SpacePointIterator end() const {
      return SpacePointIterator(
          *m_spacePointContainer,
          m_spacePointIndices.data() + m_spacePointIndices.size());
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    std::span<const SpacePointIndex2> m_spacePointIndices;
  };

  SpacePointRange spacePoints(
      const SpacePointContainer2 &spacePointContainer) const {
    return SpacePointRange(spacePointContainer,
                           m_container->spacePointIndices(m_index));
  }

 private:
  ContainerType *m_container{nullptr};
  IndexType m_index{0};
};

inline MutableSeedProxy2 SeedContainer2::createSeed(
    std::span<const SpacePointIndex2> spacePoints) {
  m_entries.emplace_back(spacePoints.size(), m_spacePoints.size(),
                         -std::numeric_limits<float>::infinity(), 0.f);
  m_spacePoints.insert(m_spacePoints.end(), spacePoints.begin(),
                       spacePoints.end());
  return MutableProxyType(*this, size() - 1);
}

inline MutableSeedProxy2 SeedContainer2::at(IndexType index) {
  return MutableProxyType(*this, index);
}

inline ConstSeedProxy2 SeedContainer2::at(IndexType index) const {
  return ConstProxyType(*this, index);
}

}  // namespace Acts::Experimental
