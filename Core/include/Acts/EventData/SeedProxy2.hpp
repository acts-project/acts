// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <cassert>

namespace Acts {

class SeedContainer2;

/// A proxy class for accessing individual seeds.
template <bool read_only>
class SeedProxy2 {
 public:
  /// Indicates whether this seed proxy is read-only or if it can be
  /// modified
  static constexpr bool ReadOnly = read_only;

  /// Type alias for seed index type
  using IndexType = SeedIndex2;

  /// Type alias for container type (const if read-only)
  using ContainerType = const_if_t<ReadOnly, SeedContainer2>;

  /// Constructs a seed proxy for the given container and index.
  /// @param container The container holding the seed.
  /// @param index The index of the seed in the container.
  SeedProxy2(ContainerType &container, IndexType index) noexcept
      : m_container{&container}, m_index{index} {}

  /// Copy construct a seed proxy.
  /// @param other The seed proxy to copy.
  SeedProxy2(const SeedProxy2 &other) noexcept = default;

  /// Copy construct a mutable seed proxy.
  /// @param other The mutable seed proxy to copy.
  explicit SeedProxy2(const SeedProxy2<false> &other) noexcept
    requires ReadOnly
      : m_container(&other.container()), m_index(other.index()) {}

  SeedProxy2 &operator=(const SeedProxy2 &other) noexcept = default;

  SeedProxy2 &operator=(const SeedProxy2<false> &other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  SeedProxy2 &operator=(SeedProxy2 &&) noexcept = default;

  SeedProxy2 &operator=(SeedProxy2<false> &&other) noexcept
    requires ReadOnly
  {
    m_container = &other.container();
    m_index = other.index();
    return *this;
  }

  /// Gets the container holding the seed.
  /// @return A reference to the container holding the seed.
  SeedContainer2 &container() noexcept
    requires(!ReadOnly)
  {
    return *m_container;
  }
  /// Gets the container holding the seed.
  /// @return A const reference to the container holding the seed.
  const SeedContainer2 &container() const noexcept { return *m_container; }
  /// Gets the index of the seed in the container.
  /// @return The index of the seed in the container.
  IndexType index() const noexcept { return m_index; }

  /// Assigns space point indices to the seed at the given index.
  /// @param spacePointIndices A span of space point indices to assign to the seed.
  /// @throws std::out_of_range if the index is out of range.
  /// @throws std::logic_error if space point indices are already assigned to the seed.
  void assignSpacePointIndices(
      std::span<const SpacePointIndex2> spacePointIndices)
    requires(!ReadOnly)
  {
    m_container->assignSpacePointIndices(m_index, spacePointIndices);
  }

  /// Returns the size of the seed, i.e., the number of space points
  /// associated with it.
  /// @return The number of space points in the seed.
  [[nodiscard]] std::size_t size() const noexcept {
    return m_container->spacePointIndices(m_index).size();
  }
  /// Checks if the seed is empty, i.e., has no space points associated with it.
  /// @return True if the seed is empty, false otherwise.
  [[nodiscard]]
  bool empty() const noexcept {
    return size() == 0;
  }

  /// Mutable access to the space point indices of the seed.
  /// @return A mutable span of space point indices associated with the seed.
  std::span<SpacePointIndex2> spacePointIndices() noexcept
    requires(!ReadOnly)
  {
    return m_container->spacePointIndices(m_index);
  }
  /// Mutable access to the quality of the seed.
  /// @return A mutable reference to the quality of the seed.
  float &quality() noexcept
    requires(!ReadOnly)
  {
    return m_container->quality(m_index);
  }
  /// Mutable access to the vertex Z coordinate of the seed.
  /// @return A mutable reference to the vertex Z coordinate of the seed.
  float &vertexZ() noexcept
    requires(!ReadOnly)
  {
    return m_container->vertexZ(m_index);
  }

  /// Const access to the space point indices of the seed.
  /// @return A span of space point indices associated with the seed.
  ///         This span is read-only and cannot be modified.
  std::span<const SpacePointIndex2> spacePointIndices() const noexcept {
    return m_container->spacePointIndices(m_index);
  }
  /// Const access to the quality of the seed.
  /// @return The quality of the seed.
  float quality() const noexcept { return m_container->quality(m_index); }
  /// Const access to the vertex Z coordinate of the seed.
  /// @return The vertex Z coordinate of the seed.
  float vertexZ() const noexcept { return m_container->vertexZ(m_index); }

  /// Iterator over space points referenced by the seed.
  class SpacePointIterator {
   public:
    /// Iterator value type
    using value_type = ConstSpacePointProxy2;
    /// Iterator difference type
    using difference_type = std::ptrdiff_t;
    /// Iterator pointer type
    using pointer = void;
    /// Iterator reference type
    using reference = void;

    /// Iterator category
    using iterator_category = std::random_access_iterator_tag;
    /// Iterator concept
    using iterator_concept = std::random_access_iterator_tag;

    SpacePointIterator() = default;
    /// Constructor from space point container and index pointer
    /// @param spacePointContainer Container of space points
    /// @param indexPointer Pointer to space point index
    SpacePointIterator(const SpacePointContainer2 &spacePointContainer,
                       const SpacePointIndex2 *indexPointer) noexcept
        : m_spacePointContainer{&spacePointContainer},
          m_indexPointer{indexPointer} {}

    /// Dereference operator
    /// @return Proxy to the space point at the current position
    value_type operator*() const noexcept {
      return (*m_spacePointContainer)[*m_indexPointer];
    }
    /// Subscript operator
    /// @param n Offset from the current position
    /// @return Proxy to the space point at offset n
    value_type operator[](difference_type n) const noexcept {
      return (*m_spacePointContainer)[m_indexPointer[n]];
    }

    /// Pre-increment operator
    /// @return Reference to the incremented iterator
    constexpr SpacePointIterator &operator++() noexcept {
      ++m_indexPointer;
      return *this;
    }
    /// Post-increment operator
    /// @return Copy of the iterator before increment
    constexpr SpacePointIterator operator++(int) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }
    /// Pre-decrement operator
    /// @return Reference to the decremented iterator
    constexpr SpacePointIterator &operator--() noexcept {
      --m_indexPointer;
      return *this;
    }
    /// Post-decrement operator
    /// @return Copy of the iterator before decrement
    constexpr SpacePointIterator operator--(int) noexcept {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    /// Compound addition assignment operator
    /// @param n Number of positions to advance
    /// @return Reference to the advanced iterator
    constexpr SpacePointIterator &operator+=(difference_type n) noexcept {
      m_indexPointer += n;
      return *this;
    }
    /// Compound subtraction assignment operator
    /// @param n Number of positions to move back
    /// @return Reference to the moved iterator
    constexpr SpacePointIterator &operator-=(difference_type n) noexcept {
      m_indexPointer -= n;
      return *this;
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    const SpacePointIndex2 *m_indexPointer{nullptr};

    friend constexpr SpacePointIterator operator+(SpacePointIterator it,
                                                  difference_type n) noexcept {
      return it += n;
    }

    friend constexpr SpacePointIterator operator+(
        difference_type n, SpacePointIterator it) noexcept {
      return it += n;
    }

    friend constexpr SpacePointIterator operator-(SpacePointIterator it,
                                                  difference_type n) noexcept {
      return it -= n;
    }

    friend constexpr difference_type operator-(
        const SpacePointIterator &lhs, const SpacePointIterator &rhs) noexcept {
      return lhs.m_indexPointer - rhs.m_indexPointer;
    }

    friend constexpr auto operator<=>(const SpacePointIterator &a,
                                      const SpacePointIterator &b) noexcept {
      return a.m_indexPointer <=> b.m_indexPointer;
    }
    friend constexpr bool operator==(const SpacePointIterator &a,
                                     const SpacePointIterator &b) noexcept {
      return a.m_indexPointer == b.m_indexPointer;
    }
  };

  /// Range facade for the seed space points.
  class SpacePointRange {
   public:
    using size_type = std::size_t;
    using value_type = ConstSpacePointProxy2;

    /// Constructor
    /// @param spacePointContainer The space point container
    /// @param spacePointIndices The space point indices
    SpacePointRange(
        const SpacePointContainer2 &spacePointContainer,
        std::span<const SpacePointIndex2> spacePointIndices) noexcept
        : m_spacePointContainer{&spacePointContainer},
          m_spacePointIndices{spacePointIndices} {}

    /// Get the number of space points in the range
    /// @return Number of space points
    std::size_t size() const noexcept { return m_spacePointIndices.size(); }
    /// Check if the range is empty
    /// @return True if the range is empty
    bool empty() const noexcept { return size() == 0; }

    /// Subscript operator
    /// @param index Index of the space point
    /// @return Proxy to the space point at the given index
    ConstSpacePointProxy2 operator[](std::size_t index) const noexcept {
      return (*m_spacePointContainer)[m_spacePointIndices[index]];
    }

    /// Get iterator to the beginning
    /// @return Iterator to the first space point
    SpacePointIterator begin() const noexcept {
      return SpacePointIterator(*m_spacePointContainer,
                                m_spacePointIndices.data());
    }
    /// Get iterator to the end
    /// @return Iterator past the last space point
    SpacePointIterator end() const noexcept {
      return SpacePointIterator(
          *m_spacePointContainer,
          m_spacePointIndices.data() + m_spacePointIndices.size());
    }

   private:
    const SpacePointContainer2 *m_spacePointContainer{nullptr};
    std::span<const SpacePointIndex2> m_spacePointIndices;
  };

  /// Get the space points associated with this seed. The space point container
  /// is taken from the seed container.
  /// @return Range of space points for this seed
  SpacePointRange spacePoints() const {
    return SpacePointRange(m_container->spacePointContainer(),
                           m_container->spacePointIndices(m_index));
  }

  /// Get the space points associated with this seed using an external space
  /// point container.
  /// @param spacePointContainer External container holding all space points
  /// @return Range of space points for this seed
  SpacePointRange spacePoints(
      const SpacePointContainer2 &spacePointContainer) const noexcept {
    return SpacePointRange(spacePointContainer,
                           m_container->spacePointIndices(m_index));
  }

 private:
  ContainerType *m_container{nullptr};
  IndexType m_index{0};
};

}  // namespace Acts
