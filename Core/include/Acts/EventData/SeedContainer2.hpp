// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/ContainerIterator.hpp"

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
  using Index = SeedIndex2;
  using MutableProxy = MutableSeedProxy2;
  using ConstProxy = ConstSeedProxy2;

  /// Constructs and empty seed container.
  SeedContainer2() noexcept;

  /// Constructs a copy of the given seed container.
  /// @param other The seed container to copy.
  SeedContainer2(const SeedContainer2 &other) noexcept;

  /// Move constructs a seed container.
  /// @param other The seed container to move.
  SeedContainer2(SeedContainer2 &&other) noexcept;

  /// Detructs the seed container.
  ~SeedContainer2() noexcept = default;

  /// Assignment operator for copying a seed container.
  /// @param other The seed container to copy.
  /// @return A reference to this seed container.
  SeedContainer2 &operator=(const SeedContainer2 &other) noexcept;

  /// Move assignment operator for a seed container.
  /// @param other The seed container to move.
  /// @return A reference to this seed container.
  SeedContainer2 &operator=(SeedContainer2 &&other) noexcept;

  /// Returns the size of the seed container, i.e., the number of seeds
  /// contained in it.
  /// @return The number of seeds in the container.
  std::size_t size() const noexcept { return m_size; }
  /// Checks if the seed container is empty.
  /// @return True if the container is empty, false otherwise.
  bool empty() const noexcept { return size() == 0; }

  /// Reserves space for the given number of seeds.
  /// @param size The number of seeds to reserve space for.
  /// @param averageSpacePoints The average number of space points per seed.
  void reserve(std::size_t size, float averageSpacePoints = 3) noexcept;

  /// Clears the seed container, removing all seeds and space points.
  void clear() noexcept;

  /// Creates a new seed.
  /// @return A mutable proxy to the newly created seed.
  MutableProxy createSeed() noexcept;

  /// Returns a mutable proxy to the seed at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the seed to access.
  /// @return A mutable proxy to the seed at the given index.
  /// @throws std::out_of_range if the index is out of range.
  MutableProxy at(Index index);
  /// Returns a const proxy to the seed at the given index.
  /// If the index is out of range, an exception is thrown.
  /// @param index The index of the seed to access.
  /// @return A const proxy to the seed at the given index.
  /// @throws std::out_of_range if the index is out of range.
  ConstProxy at(Index index) const;

  /// Returns a mutable proxy to the seed at the given index.
  /// @param index The index of the seed to access.
  /// @return A mutable proxy to the seed at the given index.
  MutableProxy operator[](Index index) noexcept;
  /// Returns a const proxy to the seed at the given index.
  /// @param index The index of the seed to access.
  /// @return A const proxy to the seed at the given index.
  ConstProxy operator[](Index index) const noexcept;

  /// Assigns space point indices to the seed at the given index.
  /// @param index The index of the seed to assign space point indices to.
  /// @param spacePointIndices A span of space point indices to assign to the seed.
  /// @throws std::out_of_range if the index is out of range.
  /// @throws std::logic_error if space point indices are already assigned to the seed.
  void assignSpacePointIndices(
      Index index, std::span<const SpacePointIndex2> spacePointIndices);

  /// Mutable access to the space point indices of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A span of space point indices associated with the seed at the given
  ///         index.
  std::span<SpacePointIndex2> spacePointIndices(Index index) noexcept {
    assert(index < m_spacePointCounts.size() && "Index out of bounds");
    assert(index < m_spacePointOffsets.size() && "Index out of bounds");
    return std::span<SpacePointIndex2>(
        m_spacePoints.data() + m_spacePointOffsets[index],
        m_spacePoints.data() + m_spacePointOffsets[index] +
            m_spacePointCounts[index]);
  }

  /// Mutable access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the quality of the seed at the given index.
  float &quality(Index index) noexcept {
    assert(index < m_qualities.size() && "Index out of bounds");
    return m_qualities[index];
  }
  /// Mutable access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed at the
  float &vertexZ(Index index) noexcept {
    assert(index < m_vertexZs.size() && "Index out of bounds");
    return m_vertexZs[index];
  }

  /// Const access to the space point indices of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A span of space point indices associated with the seed at the given
  ///         index.
  std::span<const SpacePointIndex2> spacePointIndices(
      Index index) const noexcept {
    assert(index < m_spacePointCounts.size() && "Index out of bounds");
    assert(index < m_spacePointOffsets.size() && "Index out of bounds");
    return std::span<const SpacePointIndex2>(
        m_spacePoints.data() + m_spacePointOffsets[index],
        m_spacePoints.data() + m_spacePointOffsets[index] +
            m_spacePointCounts[index]);
  }

  /// Const access to the quality of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the quality of the seed at the given index.
  float quality(Index index) const noexcept {
    assert(index < m_qualities.size() && "Index out of bounds");
    return m_qualities[index];
  }
  /// Const access to the vertex Z coordinate of the seed at the given index.
  /// @param index The index of the seed.
  /// @return A const reference to the vertex Z coordinate of the seed at the
  ///         given index.
  float vertexZ(Index index) const noexcept {
    assert(index < m_vertexZs.size() && "Index out of bounds");
    return m_vertexZs[index];
  }

  template <bool read_only>
  using Iterator = ContainerIterator<
      SeedContainer2,
      std::conditional_t<read_only, ConstSeedProxy2, MutableSeedProxy2>, Index,
      read_only>;

  using iterator = Iterator<false>;
  using const_iterator = Iterator<true>;

  iterator begin() noexcept { return iterator(*this, 0); }
  iterator end() noexcept { return iterator(*this, size()); }

  const_iterator begin() const noexcept { return const_iterator(*this, 0); }
  const_iterator end() const noexcept { return const_iterator(*this, size()); }

 private:
  std::uint32_t m_size{0};
  std::vector<std::size_t> m_spacePointOffsets;
  std::vector<std::uint8_t> m_spacePointCounts;
  std::vector<float> m_qualities;
  std::vector<float> m_vertexZs;
  std::vector<SpacePointIndex2> m_spacePoints;

  auto knownColumns() & noexcept {
    return std::tie(m_spacePointOffsets, m_spacePointCounts, m_qualities,
                    m_vertexZs, m_spacePoints);
  }
  auto knownColumns() const & noexcept {
    return std::tie(m_spacePointOffsets, m_spacePointCounts, m_qualities,
                    m_vertexZs, m_spacePoints);
  }
  auto knownColumns() && noexcept {
    return std::tuple(std::move(m_spacePointOffsets),
                      std::move(m_spacePointCounts), std::move(m_qualities),
                      std::move(m_vertexZs), std::move(m_spacePoints));
  }
};

}  // namespace Acts::Experimental

#include "Acts/EventData/SeedContainer2.ipp"
